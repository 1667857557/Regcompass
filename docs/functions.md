# Function reference

## Complete workflow

- `rc_run_regcompass_one_shot()`: prepare species-aware defaults and run the complete workflow.
- `rc_run_regcompass()`: run the complete workflow with explicit GEM, media, stage arguments, and worker counts.

## Restartable stages

- `rc_regcompass_step_grn()`: filter cells, select standard or condition-specific Pando per broad cell type, and fit GRNs.
- `rc_regcompass_step_metacells()`: construct condition-pure multimodal SuperCell metacells.
- `rc_regcompass_step_meta_modules()`: build condition-by-cell-type biological reaction catalogues and cell-type unions.
- `rc_regcompass_step_layer1()`: combine RNA and regulatory support and apply Boolean GPR rules.
- `rc_regcompass_step_layer2()`: reconstruct one cell-type-by-medium model with original MATLAB CORDA2 by default and score directional reactions. FASTCORE and the COMPASS-style complete full GEM are explicit supplementary routes. Medium scenarios modify exchange bounds only and do not directly remove reactions.
- `rc_regcompass_step_results()`: assemble annotations, evidence classes, rankings, metacell tables, and condition contrasts.

The Layer 2 parameters and alternatives are documented in [Layer 2 model builders](layer2-model-builders.md).

## GEM and medium preparation

- `rc_prepare_gem()`: load and prepare a supported human or mouse GEM.
- `rc_prepare_human2_gem()`: prepare the bundled or downloaded Human-GEM 2 model.
- `rc_prepare_mouse_gem()`: prepare the supported mouse GEM.
- `rc_make_medium_scenarios()`: construct built-in or custom medium tables.
- `rc_bundled_gem_manifest()`: inspect bundled model availability.
- `rc_download_species_gem()`: download a supported species model.

## Post analysis

- `rc_regcompass_step_target_union()`: use selected core reactions or genes as anchors, identify directly database-linked non-core reactions, and rescore them in the exact existing cell-type structural models without rebuilding Layer 2.
- `rc_test_condition_reactions()`: perform pairwise Wilcoxon and optional Kruskal-Wallis condition comparisons for fixed cell type, reaction direction, and medium.
- `rc_plot_condition_reaction()`: draw violin, violin-plus-boxplot, or boxplot distributions for one selected reaction target across conditions, with metacell points and significance annotations.
- `rc_select_gene_reactions()`: select reactions through Boolean GPR annotations for specified metabolic genes.
- `rc_plot_condition_gene_reactions()`: test and plot significant reaction directions associated with specified metabolic genes.
- `rc_build_reaction_annotations()`: create formal reaction names, formulas, substrates, products, GPRs, and database identifiers.
- `rc_attach_reaction_annotations()`: attach reaction annotations and evidence provenance to an existing result.

The optional limma metacell-level differential workflow is documented in [Post analysis](tutorial-04-post-analysis.md). It is intentionally described as within-dataset metacell inference rather than donor-level biological-replicate inference.

## Export and execution

- `rc_parallel_config()`: inspect the resolved platform-aware parallel configuration.

Use the generated Rd help for complete argument definitions. Principles and equations are maintained in [mathematical-model.md](mathematical-model.md).
