# RegCompassR

RegCompassR connects paired single-cell RNA+ATAC regulatory evidence to
metacell-level GEM and COMPASS-like reaction scoring.

## Minimal workflow

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L
  ),
  layer1_args = list(
    gpr_and_method = "min"
  ),
  medium_scenarios = medium_scenarios,
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  )
)
```

The canonical metacell default is `gamma = 30`, corresponding to an approximate
target of 30 cells per parent metacell before the post-clustering condition
split.

Built-in biological scenarios are:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

The human nutrient composition is anchored to HPLM from *Cell* 2017 and its
updated formulation from *Cell Metabolism* 2021. Plasmax from *Science
Advances* 2019 is retained as independent validation and is not numerically
averaged with HPLM. `mouse_plasma` uses a conservative metabolite set anchored
to absolute mouse plasma and interstitial-fluid measurements; unsupported
components are omitted rather than copied from human HPLM.

All five culture challenges use the same HPLM 2017/2021 basal composition and
override only the named glucose, lactate, or glutamine concentration from the
challenge paper. Background and challenge references are stored separately in
the output. User-defined reaction or metabolite compositions remain supported
through `scenario = "custom"` or `scenario = NULL` with `custom_medium` or
`custom_metabolites`.

With at least two retained conditions, Stage 1 uses `condition_grn`; otherwise
it uses `standard_pando` through `Pando::infer_grn()` and calculates no condition
coefficients.

A canonical run may explicitly omit condition metadata:

```r
single_result <- rc_run_regcompass(
  object = A,
  gem = gem,
  outdir = "RegCompass_single",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = NULL,
  celltype_col = "cell_type"
)
```

In this route, `single_result$reaction_ranking` remains available and
`single_result$condition_contrast` is empty.

## Condition-aware regulatory model

For every broad cell type with at least two retained conditions, Pando performs:

1. biological candidate discovery in the pooled cell type and separately in
   every condition;
2. exact union of observed `(TF, peak, target)` triples, without a Cartesian
   product;
3. one frozen edge dictionary shared by all conditions of that cell type;
4. one Gaussian identity GLM per condition using the same unscaled predictor
   definitions and TF–peak interaction term;
5. within-condition BH adjustment over the complete frozen dictionary.

RegCompass accepts a regulatory edge only when it is estimable and has
`padj < 0.05`. Its `penalty_effect` is the fitted condition coefficient for an
accepted edge and zero otherwise. Non-estimable coefficients remain `NA` in the
complete coefficient table and are not interpreted as biological zeros. The
pooled fit is used for candidate recall only; coefficients are not rescaled by a
pooled coefficient.

Historical output names containing `condition_full_oof`, `common`, or
OOF and does not fit a shared-slope/condition-deviation decomposition.

Stage 2 calls `SuperCell::SCimplify_by_graph_group()` with
`cell.graph.group = cell_type` and `cell.split.condition = condition`. For each
broad cell type, all conditions jointly determine one native RNA+ATAC WNN graph,
adaptive modality weights, neighbours, and Walktrap parent clusters. Condition
splits parent membership only after clustering, yielding condition-pure final
metacells without condition-specific graph fitting. No sample-derived grouping
or concatenated condition-by-cell-type stratum is used.

## Cell-type structural models

Stage 3 first constructs condition-specific biological meta-modules, then unions
them **only within the same cell type**. Different cell types retain independent
core and reaction-membership catalogues.

For `model_mode = "meta_module_gem"`, Stage 5 builds one structural model for
every `cell_type × medium_scenario` pair. FASTCORE runs independently within
each of those cell-type union GEMs. Conditions and metacells of the same cell
type reuse the corresponding model; different cell types never share a union
GEM, FASTCORE support set, model checksum, or directional `vmax` cache.

The primary metabolic ranking uses the condition-specific fixed-dictionary
primary route, the condition-unique increment is a zero compatibility matrix,
and RNA-only scoring remains an interpretation control.

## Optional targeted reaction remapping

After a completed `model_mode = "meta_module_gem"` run,
`rc_regcompass_step_target_union()` can score direct KEGG-, Reactome-, or
master-Rhea-linked non-core reactions. It reuses the exact cached union GEMs for
the corresponding cell type and medium, does not rebuild a model, and does not
rerun FASTCORE. Candidate availability is intersected across media within one
cell type, never across cell types.

## Documentation

- [Tutorial 1: one-shot workflow](docs/tutorial-01-quick-start.md)
- [Tutorial 2: stepwise workflow](docs/tutorial-02-stepwise-audit.md)
- [Tutorial 3: mathematical model](docs/tutorial-03-mathematical-model.md)
- [Tutorial 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Tutorial 5: condition comparison](docs/tutorial-05-condition-differential-analysis.md)
- [Medium scenarios and references](docs/medium-presets.md)
- [Public functions](docs/functions.md)
- [Stage schemas](docs/stage-interface-contracts.md)

Metacells are valid within-dataset statistical units. Their P values quantify
condition-associated metacell separation and are not donor-level biological
replicate inference.
