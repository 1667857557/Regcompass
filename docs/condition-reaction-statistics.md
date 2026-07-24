# Condition-associated reaction statistics

`rc_test_condition_reactions()` compares the same reaction target across
conditions within the same cell type after Layer 2 scoring. Stage 6 now attaches
formal reaction names, stoichiometric formulas, substrates, products, GPR rules,
participating genes, and RNA-versus-multiome evidence provenance to every result
row.

## Why the comparison is valid

RegCompass builds one deduplicated global union-GEM from the locally completed
meta-modules identified across all condition-by-cell-type strata. Every metacell
is scored with the same stoichiometric matrix, bounds, medium, reaction
direction, and target-flux fraction. The only unit-specific term in the Layer 2
objective is the multiome-derived reaction penalty.

For reaction target `r` and unit `j`, the tested support score is

```text
S[r,j] = -log(P[r,j] / (omega * Vmax[r]) + eps)
```

where `P[r,j]` is the minimum evidence-discordance penalty required to sustain
`omega * Vmax[r]`. Larger scores indicate stronger support. Before testing, the
function verifies that `Vmax[r]` is invariant across units under the shared GEM.

## Reaction annotations created by Stage 6

A newly generated complete result contains:

```r
result$reaction_catalog
result$reaction_evidence
```

`reaction_catalog` has one row per scored GEM reaction and includes:

- `reaction_name`;
- `subsystem` and `reaction_role`;
- `model_formula` reconstructed from the GEM stoichiometric matrix;
- `forward_substrates`, `forward_products`, and `forward_formula`;
- `reverse_substrates`, `reverse_products`, and `reverse_formula`;
- `genes`, `gpr_rule`, and `n_gpr_genes`;
- KEGG, Reactome, Rhea, and master-Rhea identifiers when available.

Metabolite names are followed by compartments, for example:

```text
L-glutamate [c] + ATP [c] -> ADP [c] + ...
```

The statistics tables additionally contain `tested_substrates`,
`tested_products`, and `tested_formula`, which follow the scored forward or
reverse LP direction rather than only the model's canonical equation order.

For results created before this feature, attach the exact GEM used in the run:

```r
result <- rc_attach_reaction_annotations(
  result,
  gem,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem"
)
```

## RNA-only versus RNA+ATAC evidence

`reaction_evidence` contains one row per
`condition × cell type × reaction`. Evidence classes are:

```text
RNA+ATAC
RNA-only
GPR/no-observed-RNA
structural/no-GPR
```

The classification is deliberately reaction-level rather than merely gene-level.
For every GPR reaction, RegCompass computes two otherwise identical reaction
capacities with the canonical GPR aggregation rules:

```text
C_RNA       = GPR(gene_support_rna)
C_RNA+ATAC  = GPR(gene_support_multiome)
```

The evidence classes mean:

- `RNA+ATAC`: `C_RNA+ATAC` differs from `C_RNA` by more than
  `evidence_tolerance` in at least one metacell of the condition–cell-type group;
- `RNA-only`: the reaction has positive RNA-derived reaction capacity, but
  ATAC integration does not change the GPR-aggregated reaction capacity in that
  group;
- `GPR/no-observed-RNA`: a GPR exists but its RNA-only reaction capacity is not
  positively supported in the group;
- `structural/no-GPR`: the reaction has no GPR and is supported structurally,
  such as many exchange, demand, sink, or artificial-support reactions.

The evidence table retains both reaction-level and gene-level provenance:

```text
evidence_resolution
has_rna_evidence
has_atac_regulatory_evidence
has_active_multiome_contribution
rna_supported_genes
atac_modifier_genes
multiome_contributing_genes
median_multiome_capacity_shift
max_abs_multiome_capacity_shift
```

`has_atac_regulatory_evidence` means that at least one GPR gene has a non-zero
accessibility-derived regulatory modifier. This is not sufficient by itself to
label the whole reaction `RNA+ATAC`: the modifier may act on a non-limiting
isozyme or subunit, or zero-preserving RNA support may remain zero.

`has_active_multiome_contribution` is stricter. It is TRUE only when the
GPR-aggregated reaction capacity changes after RNA+ATAC integration. The two
capacity-shift columns quantify the direction and magnitude of that reaction-
level change. Canonical results use `evidence_resolution = "reaction_capacity"`.

## Basic condition test

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  min_units = 5,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium",
  outdir = "RegCompass_result/07_condition_statistics"
)
```

The enriched tables can be inspected directly:

```r
condition_stats$omnibus[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "target_direction",
    "tested_formula",
    "genes",
    "gpr_rule",
    "evidence_by_condition",
    "p_adj"
  )
]

condition_stats$pairwise[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "target_direction",
    "tested_formula",
    "condition_a",
    "condition_b",
    "evidence_class_a",
    "evidence_class_b",
    "evidence_resolution_a",
    "evidence_resolution_b",
    "multiome_contributing_genes_a",
    "multiome_contributing_genes_b",
    "median_multiome_capacity_shift_a",
    "median_multiome_capacity_shift_b",
    "delta_median_score_b_minus_a",
    "rank_biserial_b_minus_a",
    "p_adj"
  )
]
```

When `outdir` is supplied, the annotated exports include:

```text
condition_reaction_pairwise.tsv.gz
condition_reaction_omnibus.tsv.gz
condition_reaction_catalog.tsv.gz
condition_reaction_evidence.tsv.gz
condition_reaction_statistics.rds
```

## Multi-condition comparison

When three or more conditions are retained, one Kruskal-Wallis omnibus test is
run for every fixed `cell type × reaction × direction × medium` target. Pairwise
Wilcoxon tests are also returned for every requested condition pair.

`condition_stats$omnibus` answers whether at least one condition differs.
`condition_stats$pairwise` provides the direction, effect size, and evidence
provenance for each contrast.

Positive `delta_median_score_b_minus_a`, Cohen's d, or rank-biserial correlation
means stronger reaction support in `condition_b`.

## Select reactions by metabolic genes

A gene query selects reactions through the Boolean GPR annotation, not by
reaction name or subsystem text.

```r
rela_metabolic_genes <- c(
  "SLC7A11",
  "GCLC",
  "GCLM",
  "GSS",
  "GSR",
  "G6PD",
  "PGD"
)

gene_reactions <- rc_select_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  match = "any",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like"
)

gene_reactions$reactions[
  ,
  c(
    "reaction_id",
    "reaction_name",
    "model_formula",
    "genes",
    "gpr_rule",
    "matched_genes"
  )
]

gene_reactions$evidence[
  ,
  c(
    "reaction_id",
    "condition",
    "cell_type",
    "evidence_class",
    "evidence_resolution",
    "rna_supported_genes",
    "atac_modifier_genes",
    "multiome_contributing_genes",
    "median_multiome_capacity_shift",
    "max_abs_multiome_capacity_shift"
  )
]
```

To retain only reactions whose GPR-aggregated reaction capacity is actively
changed by RNA+ATAC integration in at least one selected group:

```r
multiome_gene_reactions <- rc_select_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  cell_types = "stem-cell_like",
  evidence_class = "RNA+ATAC"
)
```

Selecting a gene means that it participates in the reaction GPR. It does not
mean that the gene alone is sufficient for activity. Always inspect `gpr_rule`
for multisubunit complexes and alternative isozymes.

## Plot one annotated reaction

```r
p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "epithelial_like",
  target_direction = "reverse",
  medium_scenario = "high_glucose",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  annotation_p = "p_adj"
)

print(p)
```

The default plot title uses the formal reaction name when available. The caption
contains the tested directional formula, participating genes, and evidence class
for the compared conditions.

The plotter computes multiplicity correction over the full scored reaction set
within the selected statistical scope before extracting the requested reaction.
It therefore does not treat one displayed reaction as the complete testing
family.

## Plot a significant reaction collection for selected genes

`rc_plot_condition_gene_reactions()` runs the statistics once, selects GPR
reactions containing the requested genes, filters significant targets, ranks
them by adjusted P value and effect size, and returns one annotated boxplot per
reaction direction.

```r
gene_plots <- rc_plot_condition_gene_reactions(
  result,
  genes = rela_metabolic_genes,
  cell_type = "stem-cell_like",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  comparisons = list(
    c("control_24hr", "MS177_24hr"),
    c("JQ1_24hr", "MS177_24hr")
  ),
  target_directions = c("forward", "reverse"),
  medium_scenario = "high_glucose",
  evidence_class = "RNA+ATAC",
  p_adj_max = 0.05,
  min_abs_rank_biserial = 0.30,
  max_reactions = 12,
  outdir = "RegCompass_result/07_condition_statistics/RELA_gene_plots"
)

names(gene_plots$plots)
gene_plots$selected_targets
gene_plots$pairwise_hits

print(gene_plots$plots[[1]])
```

The returned object contains:

```text
plots
selected_targets
pairwise_hits
statistics
gene_selection
```

When `outdir` is provided, each plot is saved as PDF together with the selected
target and pairwise-result tables.

## Candidate filtering

P values should be interpreted together with effect sizes and evidence source:

```r
hits <- subset(
  condition_stats$pairwise,
  p_adj < 0.05 &
    abs(rank_biserial_b_minus_a) >= 0.30 &
    abs(delta_median_score_b_minus_a) >= 0.10 &
    (evidence_class_a == "RNA+ATAC" |
       evidence_class_b == "RNA+ATAC")
)
```

Forward and reverse targets remain separate because they are separate LP
objectives. A significant score difference indicates differential support for a
reaction direction, not direct observation of net flux.

## Inference boundary

For `unit = "metacell"`, the P values quantify within-dataset separation of
condition-associated metacell score distributions. They do not make metacells
independent biological replicates. Results therefore carry:

```text
inference_level = metacell_within_dataset
descriptive_only = TRUE
biological_replicate_inference = FALSE
```

The RNA-versus-multiome evidence label describes which measured evidence enters
the reaction support calculation. It does not convert metacells into biological
replicates and does not demonstrate actual metabolic flux. Independent samples
and targeted metabolomics or isotope tracing remain necessary for population-
level and flux-level validation.
