"""Executable oracle for resendislab/corda at the pinned CORDA2 commit.

The script intentionally validates the source as written, including the
positive-coefficient minimization in build iteration 2, the unblocked opposite
reversible variable, and forward-confidence penalty coupling.
"""

from __future__ import annotations

import inspect
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from corda import CORDA
from corda import corda as corda_module

ROOT = Path(__file__).resolve().parent
ORACLE = ROOT / "corda2_python_oracle.tsv"


def add_reaction(model: Model, rid: str, stoich: dict, bounds=(0.0, 1000.0)):
    reaction = Reaction(rid)
    reaction.bounds = bounds
    reaction.add_metabolites(stoich)
    model.add_reactions([reaction])
    return reaction


def normalized_variables(worker: CORDA):
    mapping = {}
    for reaction in worker.model.reactions:
        mapping[reaction.id] = f"{reaction.id}::forward"
        mapping[reaction.reverse_id] = f"{reaction.id}::reverse"
    return mapping


def normalized_conf(worker: CORDA):
    mapping = normalized_variables(worker)
    return {mapping[name]: int(value) for name, value in worker.conf.items()}


def normalized_impossible(worker: CORDA):
    mapping = normalized_variables(worker)
    return sorted(mapping[name] for name in worker.impossible)


def normalized_redundancies(worker: CORDA):
    mapping = normalized_variables(worker)
    return {mapping[name]: int(value) for name, value in worker.redundancies.items()}


def source_contract_checks():
    source = inspect.getsource(CORDA.associated)
    build = inspect.getsource(CORDA.build)
    init = inspect.getsource(CORDA.__init__)

    assert corda_module.UPPER == 1e6
    assert corda_module.CI == 1.01
    assert "penalties[r.forward_variable] = pen" in source
    assert "penalties[r.reverse_variable] = pen" in source
    assert "conf[r.id]" in source
    assert "va.lb = max(self.tflux, va.lb)" in source
    assert "va.ub = UPPER" in source
    assert "reverse_variable" not in source.split("for vid in targets:", 1)[1]
    assert "self.model.objective.set_linear_coefficients({v: 1})" in build
    assert "self.model.objective.value > self.tflux" in build
    assert "self.model.variables[vid].ub = 0.0" in build
    assert "self.tflux = 1" in init


def association_redundancy_case(rows):
    a, b, c = Metabolite("A"), Metabolite("B"), Metabolite("C")
    model = Model("association_redundancy")
    add_reaction(model, "R1", {a: -1, c: 1})
    add_reaction(model, "R2", {b: -1, c: 1})
    add_reaction(model, "SRC_A", {a: 1})
    add_reaction(model, "SRC_B", {b: 1})
    add_reaction(model, "SINK_C", {c: -1})
    confidence = {
        "R1": 1,
        "R2": 2,
        "SRC_A": 1,
        "SRC_B": 1,
        "SINK_C": 3,
    }
    worker = CORDA(model, confidence, n=3)
    mapping = normalized_variables(worker)
    needed = worker.associated(["SINK_C"], redundancies=True)
    normalized = sorted(mapping[name] for name in needed)
    expected = sorted([
        "R1::forward", "R2::forward", "SRC_A::forward", "SRC_B::forward"
    ])
    assert normalized == expected, normalized
    assert worker.redundancies["SINK_C"] == 2
    rows.append(("association", "needed", ";".join(normalized)))
    rows.append(("association", "redundancy", "2"))


def reversible_self_cycle_case(rows):
    a, b = Metabolite("A"), Metabolite("B")
    model = Model("reversible_self_cycle")
    add_reaction(model, "REV", {a: -1, b: 1}, bounds=(-1000.0, 1000.0))
    worker = CORDA(model, {"REV": 3}, n=1)
    needed = worker.associated(["REV"], redundancies=False)
    assert len(needed) == 0
    assert "REV" not in worker.impossible
    assert worker.model.variables["REV"].primal > worker.tol
    reverse_id = worker.model.reactions.get_by_id("REV").reverse_id
    assert worker.model.variables[reverse_id].primal > worker.tol
    rows.append(("self_cycle", "status", "optimal"))
    rows.append(("self_cycle", "opposite_blocked", "false"))


def build_stage3_case(rows):
    a = Metabolite("A")
    model = Model("stage3_unknown")
    add_reaction(model, "SRC", {a: 1})
    add_reaction(model, "H", {a: -1})
    worker = CORDA(model, {"SRC": 0, "H": 3}, n=3)
    worker.build()
    conf = normalized_conf(worker)
    assert conf["SRC::forward"] == 3
    assert conf["SRC::reverse"] == -1
    assert conf["H::forward"] == 3
    assert conf["H::reverse"] == -1
    for variable, value in conf.items():
        rows.append(("stage3_unknown", variable, str(value)))
    rows.append(("stage3_unknown", "included", ";".join(
        sorted(rid for rid, included in worker.included.items() if included)
    )))
    rows.append(("stage3_unknown", "impossible", ";".join(
        normalized_impossible(worker)
    )))


def build_support_case(rows):
    x = Metabolite("X")
    model = Model("absent_support")
    add_reaction(model, "N", {x: 1})
    add_reaction(model, "M1", {x: -1}, bounds=(2.0, 1000.0))
    add_reaction(model, "M2", {x: -1})
    worker = CORDA(model, {"N": -1, "M1": 2, "M2": 2}, n=1, support=2)
    worker.build()
    conf = normalized_conf(worker)
    assert conf["N::forward"] == 3
    assert conf["N::reverse"] == -1
    assert conf["M1::forward"] == 3
    assert conf["M2::forward"] == 2
    for variable, value in conf.items():
        rows.append(("absent_support", variable, str(value)))
    rows.append(("absent_support", "included", ";".join(
        sorted(rid for rid, included in worker.included.items() if included)
    )))


def build_positive_min_case(rows, forced: bool):
    a = Metabolite("A")
    model = Model("positive_min_forced" if forced else "positive_min_free")
    add_reaction(model, "SRC", {a: 1})
    bounds = (2.0, 1000.0) if forced else (0.0, 1000.0)
    add_reaction(model, "M", {a: -1}, bounds=bounds)
    worker = CORDA(model, {"SRC": 0, "M": 2}, n=1)
    worker.build()
    conf = normalized_conf(worker)
    case = "positive_min_forced" if forced else "positive_min_free"
    if forced:
        assert conf["M::forward"] == 3
        assert conf["SRC::forward"] == 3
    else:
        assert conf["M::forward"] == 2
        assert conf["SRC::forward"] == -1
    for variable, value in conf.items():
        rows.append((case, variable, str(value)))
    rows.append((case, "included", ";".join(
        sorted(rid for rid, included in worker.included.items() if included)
    )))


def main():
    source_contract_checks()
    rows = []

    # Solver configuration is part of the source semantics.
    probe = Model("tolerance_probe")
    a = Metabolite("A")
    add_reaction(probe, "R", {a: 1})
    worker = CORDA(probe, {"R": 3})
    rows.append(("constants", "UPPER", str(corda_module.UPPER)))
    rows.append(("constants", "CI", str(corda_module.CI)))
    rows.append(("constants", "tflux", str(worker.tflux)))
    rows.append(("constants", "feasibility_tolerance", repr(worker.tol)))
    rows.append(("constants", "default_n", str(worker.n)))
    rows.append(("constants", "default_penalty_factor", str(worker.pf)))
    rows.append(("constants", "default_support", str(worker.support)))

    association_redundancy_case(rows)
    reversible_self_cycle_case(rows)
    build_stage3_case(rows)
    build_support_case(rows)
    build_positive_min_case(rows, forced=True)
    build_positive_min_case(rows, forced=False)

    with ORACLE.open("w", encoding="utf-8") as handle:
        handle.write("case\tkey\tvalue\n")
        for case, key, value in rows:
            handle.write(f"{case}\t{key}\t{value}\n")

    print("Pinned Python CORDA2 exact-source oracle passed")


if __name__ == "__main__":
    main()
