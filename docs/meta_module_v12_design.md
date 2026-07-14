# RegCompassR 1.2: sample-specific Pando GRN meta-modules

## 1. Scope

RegCompassR 1.2 adds a sample-specific regulatory layer between metacell construction and reaction-set/GEM analysis. The existing Layer 1 reaction-capacity and ATAC-confidence calculations remain unchanged. Pando results define biologically motivated reaction collections; they do not replace GPR capacity, medium constraints, mass balance, or flux optimization.

Implementation boundary: the code currently defines the biological reaction envelope **B** as GRN → metabolic gene → GPR reaction → subsystem/cross-reference-linked subsystem expansion. The intended complete local-GEM workflow continues after **B** with medium-specific full GEM construction, FASTCC consistency, core-direction precheck, add-only FASTCORE support completion, capacity-retention validation, and sample/module-aware microCOMPASS. FASTCC/FASTCORE/capacity-retention are not yet implemented in `rc_build_meta_module_gem()`; the current function uses deterministic Human-GEM support inclusion and labels those reactions separately from biological meta-module reactions.

The implemented biological-envelope stages and intended local-GEM support-completion stages are separated as follows:

```text
single-cell RNA+ATAC
→ strict condition × sample × cell-type filtering
→ RNA+ATAC metacells and metacell fragments
→ existing metacell LinkPeaks and Layer 1 capacity/confidence
→ split retained metacells by sample
→ sample-specific Pando GRN
→ significant metabolic target genes + tightly connected metabolic genes
→ GPR mapping to core reactions
→ subsystem/database/master-Rhea expansion
→ biological reaction envelope B
→ medium-specific full GEM
→ FASTCC consistent parent
→ core direction precheck
→ add-only FASTCORE support completion
→ capacity-retention validation
→ completed local GEM = B* ∪ support
→ sample/module-aware microCOMPASS
```

## 2. Version pinning

The package pins the Seurat 4-compatible fork:

```r
remotes::install_github(
  "1667857557/Pando_regcompass"
)
```

The package name remains `Pando`. Runtime validation requires the installed package metadata to point to `1667857557/Pando_regcompass`; it records but does not restrict the package version or commit SHA.

## 3. Human-GEM annotation preparation

Use `rc_prepare_human2_gem_v12()` rather than the legacy preparation helper when constructing meta-modules.

The v1.2 object retains these normalized reaction metadata columns:

| Column | Human-GEM source | Use |
|---|---|---|
| `subsystem` | reaction YAML subsystem field or supplied reaction–subsystem table | Step 3 and subsystem-level expansion |
| `kegg_reaction_id` | `rxnKEGGID` | Step 4 |
| `reactome_reaction_id` | `rxnREACTOMEID` | Step 4 |
| `rhea_reaction_id` | `rxnRheaID` | provenance |
| `rhea_master_id` | `rxnRheaMasterID` | Step 5 |

A reaction can have multiple subsystem or database annotations. Delimited values are normalized to long maps by `rc_reaction_crossref_maps()`.

`rxnKEGGID` and `rxnREACTOMEID` are reaction/event cross-references, not pathway identifiers. Step 4 therefore means: identify another subsystem containing at least one reaction with a KEGG/Reactome reaction identifier shared with the current subsystem collection, then add all reactions from that subsystem. A future pathway-level implementation would require an explicit subsystem-to-pathway crosswalk or pathway membership table. In other words, the present implementation should be described as subsystem/cross-reference expansion, not as a completed GRN–gene–reaction–subsystem–pathway expansion.

## 4. Per-sample Pando input

### 4.1 Analysis unit

Pando is run once per `sample_id`, after retaining only metacells that passed the formal pre- and post-metacell filters. Metacells remain the observations used for regression. Each sample object can contain multiple retained cell types; cell type is not pooled across samples.

### 4.2 RNA target genes

The Pando target set is:

```text
rownames(original single-cell RNA assay)
∩ Human-GEM GPR metabolic genes
∩ rownames(retained metacell RNA assay)
```

Matching is case-insensitive, while the exact metacell RNA feature names are passed to Pando.

This intentionally avoids Pando's default variable-feature-only target selection. All observed Human-GEM metabolic genes are eligible target genes.

### 4.3 ATAC and motif preparation

For each sample:

```r
obj <- subset(metacell_object, cells = sample_metacells)
obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
obj <- Signac::RunTFIDF(obj, assay = "ATAC")
grn <- Pando::initiate_grn(
  obj,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  exclude_exons = TRUE
)
grn <- Pando::find_motifs(grn, pfm = pfm, genome = genome)
grn <- Pando::infer_grn(
  grn,
  genes = metabolic_target_genes,
  method = "glm",
  tf_cor = 0.1,
  peak_cor = 0,
  adjust_method = "fdr",
  parallel = FALSE
)
```

`glm` is the default because the requested significant-edge filter requires p-values and adjusted p-values. Other Pando fitters can be passed through `pando_infer_args`, but `require_padj = TRUE` will reject output without `padj`.

### 4.4 Significant TF–peak–gene edges

`rc_extract_pando_tf_peak_gene()` exports the complete coefficient table and applies:

```text
padj ≤ 0.05
|estimate| ≥ min_abs_estimate
model rsq ≥ 0.1, when rsq is available
```

The complete and filtered tables retain at minimum:

```text
sample_id, tf, target, region, estimate, pval, padj, corr, rsq
```

## 5. Defining tightly connected metabolic genes

Pando directly estimates TF–peak–target relationships, not generic metabolic-gene-to-metabolic-gene edges. RegCompassR therefore creates an explicit metabolic projection with two edge types.

### 5.1 Shared-TF projection

For each sample, two significant metabolic target genes are connected when they share significant upstream TFs. For genes `a` and `b`:

```text
shared_tf_count(a,b) = number of significant TFs shared by a and b
projection_weight(a,b) = sum over shared TFs of min(|βTF,a|, |βTF,b|)
tf_jaccard(a,b) = |TFa ∩ TFb| / |TFa ∪ TFb|
```

Defaults:

```text
min_shared_tfs = 1
tf_jaccard ≥ 0
top_k_neighbors = 5 per gene
max_targets_per_tf = 200
```

The target cap prevents highly promiscuous TFs from creating quadratic, near-complete networks.

### 5.2 Direct metabolic-TF edges

When a significant TF is itself a Human-GEM metabolic gene, the metabolic TF and its metabolic target are connected directly. Direct edges are retained independently of shared-TF thresholds.

### 5.3 GRN components

Connected components of the retained metabolic projection become sample-specific modules:

```text
<sample_id>::GRN0001
<sample_id>::GRN0002
...
```

Node roles are recorded as:

- `significant_target`
- `metabolic_tf_neighbor`
- `target_and_metabolic_tf`

This definition is deterministic for a fixed Pando coefficient table and parameter set.

## 6. Mapping GRN genes to core reactions

`rc_map_meta_module_core_reactions()` joins every GRN node gene to every Human-GEM GPR reaction containing that gene.

The mapping is inclusive:

- a gene can map to multiple reactions;
- a reaction can be supported by multiple GRN genes;
- reactions containing an AND GPR are still designated core when any mapped gene places the reaction in the GRN module;
- actual reaction capacity continues to be calculated by the existing AND/OR-aware Layer 1 code.

The distinction prevents network membership from being confused with complete enzyme-complex sufficiency.

## 7. Ordered reaction expansion

`rc_expand_meta_module_reactions()` applies the requested rules independently to every `sample_id × GRN module`.

### Stage 1 — core

Include all reactions mapped from GRN metabolic genes.

```text
inclusion_stage = core_grn_gene
```

### Stage 2 — core subsystems

For every subsystem assigned to any core reaction, include all reactions assigned to that subsystem. Multi-subsystem reactions contribute all of their valid subsystem labels.

`UNASSIGNED`, `NA`, and `NONE` are never treated as shared subsystem labels.

```text
inclusion_stage = same_core_subsystem
```

### Stage 3 — KEGG/Reactome-linked subsystems

1. Collect all KEGG and Reactome reaction identifiers present among reactions currently included after Stage 2.
2. Identify every subsystem containing at least one reaction with one of these identifiers.
3. Include every reaction from each identified subsystem.

Thus, a database-linked subsystem contributes its complete reaction set, not only the cross-referenced reaction.

```text
inclusion_stage = shared_kegg_or_reactome_subsystem
```

### Stage 4 — master-Rhea-linked subsystems

1. Collect all `rhea_master_id` values among reactions currently included after Stage 3.
2. Identify reactions sharing one of these master-Rhea identifiers.
3. Identify the subsystem(s) containing those reactions.
4. Include every reaction from those subsystem(s).

```text
inclusion_stage = shared_master_rhea_subsystem
```

### Expansion mode

The default is:

```r
expansion_mode = "ordered_once"
```

This executes the four stages once and prevents uncontrolled transitive expansion.

For sensitivity analysis:

```r
expansion_mode = "fixed_point"
```

The database and Rhea rules are reapplied until no new reactions are added or `max_iterations` is reached. Fixed-point modules can become much larger and should not be treated as the primary definition without reporting both module size and expansion depth.

## 8. Biological membership versus solver support

`reaction_membership` is the exact biological reaction envelope **B** produced by GRN/gene/reaction/subsystem/cross-reference expansion. It is never silently expanded for flux feasibility.

The intended support-completion layer is:

```text
B
→ apply the requested medium to the full Human-GEM
→ run FASTCC to obtain a medium-specific consistent parent model
→ precheck each core reaction direction for achievable steady-state flux
→ run add-only FASTCORE with B/core reactions protected
→ validate capacity retention for protected reactions
→ completed local GEM = B* ∪ support
```

Here **B*** denotes the protected biological envelope after consistency/direction checks, and `support` denotes reactions selected only to make protected reactions flux-consistent. FASTCORE must be add-only: it may not delete biological reactions from **B** or promote support reactions into biological membership.

The current implementation does not yet run FASTCC, FASTCORE, or capacity-retention validation. `rc_build_meta_module_gem()` currently reuses `rc_build_module_meso_gem()` to optionally add deterministic support reactions from the complete Human-GEM:

- one-hop reactions;
- transport reactions;
- exchange reactions;
- demand/sink/maintenance/cofactor support reactions.

Every reaction in the resulting local GEM is labeled:

```text
biological_meta_module_member = TRUE/FALSE
support_only = TRUE/FALSE
```

Therefore, support reactions can satisfy local mass-balance and boundary requirements without being reported as Pando/GRN-derived module members. In the intended FASTCORE-based completion, FASTCORE-selected reactions should be labeled `support_only = TRUE` unless they were already in `reaction_membership`, and capacity-retention diagnostics should be reported before interpreting microCOMPASS scores.

Recommended primary settings:

```r
include_one_hop = FALSE
include_transport = TRUE
include_exchange = TRUE
include_protected = TRUE
```

Closure and target feasibility should still be checked before interpreting microCOMPASS scores.

## 9. Function-by-function workflow and tutorial

This section mirrors the actual implementation order. Use the integrated wrapper in Section 10 for routine runs; use the staged calls below when auditing, debugging, or rerunning a single stage.

### 9.1 Prepare a v1.2 Human-GEM object

```r
gem <- rc_prepare_human2_gem_v12(version = "2.0.0")
```

`rc_prepare_human2_gem_v12()` refuses `version = "latest"` unless `allow_latest = TRUE`, removes an existing cache only when `force_download = TRUE`, downloads the requested Human-GEM reference, converts `model/Human-GEM.yml`, enriches metadata, annotates reaction roles, validates the object, saves it to `save_rds`, and finally reads it back through `rc_read_gem()`. The returned object contains the standard stoichiometric model plus `gpr_table`, `metabolic_genes`, `reaction_rules`, `genes`, `reactions`, and v1.2 model metadata.

Supporting helpers:

- `rc_enrich_humangem_v12_metadata()` first applies the legacy Human-GEM metadata enrichment, then fills missing `subsystem` and `equation` values from the YAML file when available, copies `rxnKEGGID`, `rxnREACTOMEID`, `rxnRheaID`, and `rxnRheaMasterID` from `reactions.tsv`, and mirrors subsystem labels into `metabolic_module`.
- `rc_reaction_crossref_maps()` turns semicolon- or otherwise delimited annotations into long maps for subsystem, KEGG, Reactome, Rhea, and master-Rhea lookups. It drops `UNASSIGNED`, `NA`, and `NONE` subsystem labels before expansion.

### 9.2 Run or reuse the strict Layer 1 metacell workflow

```r
layer1 <- rc_run_regcompass_multiome_metacell(
  object = object,
  gpr_table = gem$gpr_table,
  outdir = "RegCompassR_v1.2_run",
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  gamma = 150,
  adaptive_gamma = TRUE,
  min_cells_pre_metacell = 100,
  min_metacell_size = 10,
  min_metacells_post_metacell = 10,
  future_plan = "sequential"
)
```

This is the existing formal workflow. It enforces the same `condition × sample × cell type` stratum for pre-filtering, metacell construction, post-filtering, fragment aggregation, LinkPeaks, and Layer 1 capacity/confidence. The v1.2 wrapper does not reimplement this logic; it delegates to this function.

To resume the Pando stage from a completed run:

```r
mc_object <- rc_load_metacell_object_from_run(
  run_dir = "RegCompassR_v1.2_run",
  retained_metacell_ids = layer1$metacell_meta$metacell_id,
  rna_assay = "RNA",
  atac_assay = "ATAC"
)
```

`rc_load_metacell_object_from_run()` searches `01_metacells/` recursively for `metacell_object.rds`, verifies that saved metacell IDs are unique, optionally checks/subsets retained IDs, removes stored Signac fragment pointers, merges objects if more than one file is present, and returns a retained metacell Seurat object.

### 9.3 Run the sample-specific Pando meta-module stage

```r
grns <- rc_run_pando_meta_modules(
  metacell_object = mc_object,
  gem = gem,
  outdir = file.path("RegCompassR_v1.2_run", "04_pando_meta_modules"),
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sample_col = "sample_id",
  single_cell_genes = rownames(rc_get_assay_counts(object, "RNA")),
  rna_assay = "RNA",
  atac_assay = "ATAC",
  min_metacells = 20,
  pando_infer_args = list(method = "glm", tf_cor = 0.1, peak_cor = 0,
                          adjust_method = "fdr", parallel = FALSE),
  padj_threshold = 0.05,
  min_model_rsq = 0.1,
  top_k_neighbors = 5,
  min_shared_tfs = 1,
  expansion_mode = "ordered_once"
)
```

Actual behavior by internal function:

1. `rc_run_pando_meta_modules()` validates that `metacell_object` is Seurat, validates the pinned `1667857557/Pando_regcompass` installation, requires Seurat and Signac, creates output folders, and constructs the Pando target set from original single-cell genes, retained metacell RNA genes, and Human-GEM metabolic genes. It iterates over `sample_col` values; samples with fewer than `min_metacells` are recorded as `skipped_too_few_metacells`, and sample-level errors are either recorded or stopped according to `on_sample_error`.
2. For each runnable sample, the function subsets metacells, normalizes RNA, runs ATAC TF-IDF, optionally saves the sample object, calls `Pando::initiate_grn()`, `Pando::find_motifs()`, and `Pando::infer_grn()`, then optionally saves the Pando object.
3. `rc_extract_pando_tf_peak_gene()` converts `stats::coef(grn_object)` to a data frame, requires `tf`, `target`, and `region`, merges non-conflicting `Pando::gof()` columns by `target` when available, uppercases `tf` and `target`, and filters significant rows by finite absolute `estimate`, adjusted p-value, and `rsq` when present. If `padj` is absent and `require_padj = TRUE`, it stops with guidance to use a p-value-producing Pando model such as GLM.
4. `rc_project_metabolic_grn()` keeps significant rows whose targets are metabolic genes, aggregates TF-target strengths, forms shared-TF target pairs with capped targets per TF, adds direct edges when a significant TF is itself metabolic, applies `min_shared_tfs`, `min_tf_jaccard`, and per-gene `top_k` filtering, then labels connected components as `<sample_id>::GRN####`.
5. `rc_map_meta_module_core_reactions()` uppercase-joins projected GRN genes to `gem$gpr_table` and emits one row per gene–reaction match with `is_core = TRUE` and `inclusion_stage = "core_grn_gene"`. This is membership evidence, not a replacement for AND/OR-aware GPR capacity.
6. `rc_expand_meta_module_reactions()` runs independently for each sample/module pair. It filters core reactions to valid GEM reactions, includes reactions from core subsystems, includes whole subsystems linked by KEGG or Reactome reaction IDs found in the current collection, includes whole subsystems linked by shared master-Rhea IDs, and either stops after one ordered pass or repeats in `fixed_point` mode. It does not consume a separate pathway table and should not be documented as an explicit pathway expansion.
7. The stage writes status, coefficient, projection, core, membership, and summary tables, saves `pando_meta_modules.rds`, and returns the same data as a list. It stops if no significant Pando TF–peak–gene edges are available across all samples.

### 9.4 Build a local meta-module GEM

```r
module_gem <- rc_build_meta_module_gem(
  gem = gem,
  reaction_membership = grns$reaction_membership,
  sample_id = grns$meta_module_summary$sample_id[[1]],
  module_id = grns$meta_module_summary$module_id[[1]],
  medium_table = medium,
  include_one_hop = FALSE,
  include_transport = TRUE,
  include_exchange = TRUE,
  include_protected = TRUE
)
```

`rc_build_meta_module_gem()` currently selects the requested biological reaction set **B**, writes it into `reaction_meta$grn_meta_module`, delegates deterministic support expansion and closure handling to `rc_build_module_meso_gem()`, then adds `sample_id`, `grn_module_id`, `biological_meta_module_member`, and `support_only` flags. Biological membership therefore remains auditable even when transport, exchange, protected, or other support reactions are added for feasibility. The intended replacement/extension of this support layer is medium-specific full GEM → FASTCC consistent parent → core-direction precheck → add-only FASTCORE support completion → capacity-retention validation.

### 9.5 Use the integrated convenience wrapper

`rc_run_regcompass_v12()` combines Sections 9.2 and 9.3. It requires `gem$gpr_table`, forwards `metacell_args` to `rc_run_regcompass_multiome_metacell()`, reloads retained metacells from disk, forwards `pando_args` to `rc_run_pando_meta_modules()`, stores the returned list as `layer1$grn_meta_modules`, sets `schema_version = "regcompass_v1.2"`, saves `regcompass_v1.2_result.rds`, and returns the integrated object.

## 10. Integrated example

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

gem <- rc_prepare_human2_gem_v12(version = "2.0.0")

result <- rc_run_regcompass_v12(
  object = object,
  gem = gem,
  outdir = "RegCompassR_v1.2_run",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  metacell_args = list(
    gamma = 150,
    adaptive_gamma = TRUE,
    min_cells_pre_metacell = 100,
    min_metacell_size = 10,
    min_metacells_post_metacell = 10,
    future_plan = "sequential"
  ),
  pando_args = list(
    min_metacells = 20,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr",
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1,
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    expansion_mode = "ordered_once"
  )
)

meta_modules <- result$grn_meta_modules$reaction_membership
module_summary <- result$grn_meta_modules$meta_module_summary
```

To build one local GEM:

```r
module_gem <- rc_build_meta_module_gem(
  gem = gem,
  reaction_membership = meta_modules,
  sample_id = module_summary$sample_id[[1]],
  module_id = module_summary$module_id[[1]],
  medium_table = medium,
  include_one_hop = FALSE,
  include_transport = TRUE,
  include_exchange = TRUE,
  include_protected = TRUE
)
```

## 11. Output layout

```text
RegCompassR_v1.2_run/
├── 00_stratum_qc/
├── 01_metacells/
├── 02_metacell_fragments/
├── 03_linkpeaks/
├── 04_pando_meta_modules/
│   ├── pando_sample_status.tsv.gz
│   ├── pando_tf_peak_gene_all.tsv.gz
│   ├── pando_tf_peak_gene_significant.tsv.gz
│   ├── metabolic_gene_nodes.tsv.gz
│   ├── metabolic_gene_edges.tsv.gz
│   ├── core_gene_reaction.tsv.gz
│   ├── meta_module_reactions.tsv.gz
│   ├── meta_module_summary.tsv.gz
│   ├── pando_meta_modules.rds
│   ├── sample_metacell_objects/<sample>.rds
│   └── pando_objects/<sample>.rds
└── regcompass_v1.2_result.rds
```

`meta_module_reactions.tsv.gz` is the primary audit table. Each row records sample, GRN module, reaction, core status, inclusion stage, source subsystem/database identifier, and expansion mode.

## 12. Required validation before biological interpretation

At minimum report:

1. Pando status and number of retained metacells for every sample.
2. Number of tested and significant TF–peak–gene coefficients per sample.
3. GRN node/edge/module counts and sensitivity to `padj`, `rsq`, and `top_k`.
4. Core-reaction count before subsystem expansion.
5. Reaction count added at each expansion stage.
6. Ordered-once versus fixed-point module sizes.
7. Fraction of local-GEM reactions labeled `support_only`.
8. Closure, target feasibility, and medium diagnostics.
9. Sample reproducibility or consensus frequency of genes/reactions across biological replicates.

A sample-specific edge or module is not automatically a condition-level effect. Differential conclusions must be based on replicate-level recurrence or downstream statistical testing, not on pooling all metacells as independent biological replicates.
