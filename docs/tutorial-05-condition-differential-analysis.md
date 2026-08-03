# Tutorial 5: condition-level reaction results

Condition comparisons use metacells within the same cell type, reaction, direction, and medium.

## Inspect available groups

```r
unique(result$reaction_comparison_by_metacell[, c(
  "cell_type", "condition", "medium", "direction"
)])
```

Cell types routed to standard Pando because only one condition was retained contribute reaction rankings but cannot produce a within-cell-type condition contrast.

## Reaction ranking

```r
ranking <- result$reaction_ranking
ranking <- ranking[
  ranking$cell_type == "T_cell" &
  ranking$medium == "normal_human_plasma",
]
```

## Pairwise condition contrast

```r
contrast <- result$condition_contrast
contrast <- contrast[
  contrast$cell_type == "T_cell" &
  contrast$condition_1 == "Control" &
  contrast$condition_2 == "Treatment",
]
```

Review the number of metacells per condition and the feasible/evaluated flags before interpreting a comparison.

## RNA-only control

```r
result$rna_only_control_summary
result$rna_only_control_contrast
```

The RNA-only route uses the same structural models and media and is intended as an interpretation control.

## Export

```r
write.table(
  ranking,
  file = "reaction_ranking.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

Reported statistical tests describe within-dataset metacell separation. Donor-level inference requires donor-aware biological replication outside this workflow.
