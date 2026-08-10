# CORDA2 source audit pointer

The canonical CORDA2 source used for implementation review is `schultzdre/Constraint-Based-Modeling/CORDA2.m`.

RegCompass keeps the source-derived reconstruction semantics separate from its evidence adapter and downstream COMPASS-like scoring. The current implementation-specific contract, including directional handling, immutable core reactions, and removal of the historical post-reconstruction closure pass, is maintained in [mathematical-model.md](mathematical-model.md).

This file intentionally contains no duplicate algorithm derivation. Use the git history when an implementation-change audit requires the older source-to-R comparison notes.
