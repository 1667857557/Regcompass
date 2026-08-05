"""Independent executable audit of pinned resendislab CORDA source.

This file does not import or reuse RegCompass audit fixtures. It derives the
constructor signature and deterministic behavioral cases directly from the
installed Python source at commit c02e06d50606bf93f23d8f2e6d6ade0e996ca70e.
"""

from __future__ import annotations

import inspect
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from corda import CORDA
from corda import corda as corda_module

OUT = Path(__file__).resolve().parent / "corda2_independent_oracle.tsv"


def add_reaction(model, rid, stoich=None, bounds=(0.0, 1000.0)):
    reaction = Reaction(rid)
    reaction.bounds = bounds
    if stoich:
        reaction.add_metabolites(stoich)
    model.add_reactions([reaction])
    return reaction


def normalized_variable_map(worker):
    result = {}
    for reaction in worker.model.reactions:
        result[reaction.id] = f"{reaction.id}::forward"
        result[reaction.reverse_id] = f"{reaction.id}::reverse"
    return result


def constructor_signature(rows):
    signature = inspect.signature(CORDA.__init__)
    parameters = list(signature.parameters.values())
    observed_names = [parameter.name for parameter in parameters]
    expected_names = [
        "self", "model", "confidence", "met_prod", "n",
        "penalty_factor", "support"
    ]
    assert observed_names == expected_names, observed_names
    defaults = {
        "met_prod": None,
        "n": 3,
        "penalty_factor": 100,
        "support": 5,
    }
    for name, expected in defaults.items():
        assert signature.parameters[name].default == expected
        rows.append(("signature", name, repr(expected)))
    source = inspect.getsource(CORDA.__init__)
    assert "time_limit" not in signature.parameters
    assert "time_limit" not in source
    rows.append(("signature", "parameter_order", ";".join(expected_names[1:])))
    rows.append(("signature", "has_time_limit", "false"))


def bound_case(rows, name, bounds):
    model = Model(name)
    add_reaction(model, "R", bounds=bounds)
    try:
        worker = CORDA(model, {"R": 3})
    except Exception as error:  # exact constructor rejection is the contract
        rows.append((name, "status", "error"))
        rows.append((name, "error_type", type(error).__name__))
        return

    reaction = worker.model.reactions.get_by_id("R")
    rows.extend([
        (name, "status", "ok"),
        (name, "reaction_lb", repr(float(reaction.lower_bound))),
        (name, "reaction_ub", repr(float(reaction.upper_bound))),
        (name, "forward_lb", repr(float(reaction.forward_variable.lb))),
        (name, "forward_ub", repr(float(reaction.forward_variable.ub))),
        (name, "reverse_lb", repr(float(reaction.reverse_variable.lb))),
        (name, "reverse_ub", repr(float(reaction.reverse_variable.ub))),
        (name, "tolerance", repr(float(worker.tol))),
    ])


def penalty_route_case(rows, penalty_factor):
    a = Metabolite("A")
    b = Metabolite("B")
    c = Metabolite("C")
    d = Metabolite("D")
    model = Model(f"penalty_route_{penalty_factor:g}")
    add_reaction(model, "SRC", {a: 1})
    add_reaction(model, "N", {a: -1, d: 1})
    add_reaction(model, "M1", {a: -1, b: 1})
    add_reaction(model, "M2", {b: -1, c: 1})
    add_reaction(model, "M3", {c: -1, d: 1})
    add_reaction(model, "SINK", {d: -1})
    confidence = {
        "SRC": 0,
        "N": -1,
        "M1": 1,
        "M2": 1,
        "M3": 1,
        "SINK": 3,
    }
    worker = CORDA(
        model, confidence, n=1, penalty_factor=penalty_factor, support=5
    )
    mapping = normalized_variable_map(worker)
    needed = worker.associated(["SINK"], redundancies=True)
    normalized = sorted(mapping[name] for name in needed)
    key = f"penalty_factor_{penalty_factor:g}"
    rows.append((key, "needed", ";".join(normalized)))
    rows.append((key, "redundancy", str(worker.redundancies["SINK"])))


def zero_redundancy_case(rows):
    a = Metabolite("A")
    model = Model("n_zero")
    add_reaction(model, "SRC", {a: 1})
    add_reaction(model, "SINK", {a: -1})
    worker = CORDA(model, {"SRC": 1, "SINK": 3}, n=0)
    needed = worker.associated(["SINK"], redundancies=True)
    rows.append(("n_zero", "needed", ";".join(needed.tolist())))
    rows.append(("n_zero", "redundancy", str(worker.redundancies["SINK"])))


def main():
    rows = []
    assert corda_module.UPPER == 1e6
    assert corda_module.CI == 1.01
    constructor_signature(rows)

    bound_case(rows, "bounds_reversible", (-1000.0, 1000.0))
    bound_case(rows, "bounds_positive", (2.0, 1000.0))
    bound_case(rows, "bounds_negative", (-1000.0, -2.0))
    bound_case(rows, "bounds_tiny", (-5e-8, 5e-8))
    bound_case(rows, "bounds_transient_lower_error", (-2e6, -1.5e6))
    bound_case(rows, "bounds_transient_upper_error", (1.5e6, 2e6))

    penalty_route_case(rows, 0.5)
    penalty_route_case(rows, 100.0)
    zero_redundancy_case(rows)

    with OUT.open("w", encoding="utf-8") as handle:
        handle.write("case\tkey\tvalue\n")
        for case, key, value in rows:
            handle.write(f"{case}\t{key}\t{value}\n")

    print("Independent pinned Python CORDA source audit passed")


if __name__ == "__main__":
    main()
