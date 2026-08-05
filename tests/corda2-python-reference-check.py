"""Executable behavior check for the pinned resendislab/corda reference."""

from cobra import Metabolite, Model, Reaction
from corda import CORDA


def build_model(equal_medium: bool = False):
    a = Metabolite("A")
    b = Metabolite("B")
    c = Metabolite("C")

    src_a = Reaction("SRC_A")
    src_a.add_metabolites({a: 1})
    src_b = Reaction("SRC_B")
    src_b.add_metabolites({b: 1})
    r1 = Reaction("R1")
    r1.add_metabolites({a: -1, c: 1})
    r2 = Reaction("R2")
    r2.add_metabolites({b: -1, c: 1})
    sink = Reaction("SINK_C")
    sink.add_metabolites({c: -1})

    model = Model("corda2-reference")
    model.add_reactions([src_a, src_b, r1, r2, sink])
    confidence = {
        "SRC_A": 0,
        "SRC_B": 0,
        "R1": 1,
        "R2": 1 if equal_medium else -1,
        "SINK_C": 3,
    }
    return model, confidence


def main():
    model, confidence = build_model(equal_medium=False)
    worker = CORDA(model, confidence, n=1)
    assert worker.n == 1
    assert worker.support == 5
    assert worker.pf == 100
    assert worker.tflux == 1
    needed = set(worker.associated(["SINK_C"], redundancies=False))
    assert needed == {"R1"}, needed

    model, confidence = build_model(equal_medium=True)
    worker = CORDA(model, confidence, n=3)
    needed = set(worker.associated(["SINK_C"], redundancies=True))
    assert needed == {"R1", "R2"}, needed
    assert worker.redundancies["SINK_C"] == 2
    print("Pinned Python CORDA2 reference checks passed")


if __name__ == "__main__":
    main()
