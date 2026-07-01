# RegCompassR-SeuratV4：基于 Biomni 批判报告核验后的修正版开发方案

> 输入对象：已经完成细胞分群/注释的 **Seurat v4 single-cell multiome 对象**，包含同一细胞/同一核的 RNA + ATAC 数据。  
> 开发环境：Linux + R。耗时步骤通过 `BiocParallel` 并行。  
> 工具定位：**multiome-supported, GPR-aware, sample-aware reaction capacity potential framework**。  
> 第一阶段目标：先实现稳定、可诊断、可复现的 **Layer 1 reaction capacity potential**，暂不直接进入 Human-GEM 全量 QP/FVA。

---

## 0. 对 Biomni 批判报告的核验结论

Biomni 批判报告指出的问题多数成立，尤其是以下问题必须吸收到开发方案中：

| 问题 | 核验判断 | 本方案调整 |
|---|---|---|
| `local_state` 生物学含义不清 | 成立 | `state_col` 必须显式指定来源；默认不强制使用；所有结果报告 state 来源与敏感性 |
| pool 随机划分未纳入不确定性 | 成立 | 增加 `pool_seed_replicates`，至少用于 Layer 1 sensitivity |
| pool-level TF-IDF 会低估广泛开放的 housekeeping/metabolic peaks | 基本成立 | ATAC confidence 默认改用 pool-level normalized accessibility，不用 TF-IDF 作为代谢基因调控支持主输入 |
| dropout correction 存在循环依赖与尺度不匹配 | 成立 | 主分析删除 dropout correction 对 gene score 的直接修正；detection rate 仅进入 confidence / low-power flag |
| safe MAD 量纲混用 | 成立 | 使用 `MAD(constant=1.4826)` 与 `IQR/1.349`，二者均作为 σ 估计 |
| τ=0.08 实际接近 hard min | 成立 | 不再把 τ=0.08 描述为一般“soft AND”；作为 strict-bottleneck sensitivity，主分析提供 τ grid |
| Q95 在 n<20 处仍有硬切换 | 成立 | 所有 n 使用连续 shrinkage；低 pool 数只打 `low_power_q95` flag |
| gene confidence 可能为负 | 成立 | confidence 重构为非负分量；负 RNA-ATAC 相关单独标记为 discordance |
| Fisher shrinkage 对 Spearman 的适用性有限 | 成立 | 标记为近似；n pool 太小时不给强解释 |
| Layer 1 capacity 作为 hard flux bound 生物学假设过强 | 成立 | QP 默认只使用 soft penalty；hard bound 仅作为 sensitivity |

因此，当前开发策略调整为：

```text
Seurat v4 annotated multiome
→ raw-count sample-aware pseudobulk
→ pool-level RNA/ATAC normalization
→ GPR Layer 1 capacity potential
→ multiome confidence diagnostics
→ sample-aware summary
```

第一版不做：

```text
full Human-GEM QP
FVA
thermodynamic constraints
causal regulator discovery
真实 flux 推断
```

---

## 1. 项目命名与目录结构

建议作为 R package 开发：

```text
RegCompassR/
  DESCRIPTION
  NAMESPACE
  R/
    input_seurat_v4.R
    metadata_check.R
    pooling.R
    pseudobulk.R
    normalize_pool.R
    gene_score.R
    gpr_parser.R
    gpr_capacity.R
    calibration_q95.R
    multiome_confidence.R
    diagnostics.R
    parallel.R
    sample_summary.R
    report.R
    qp_placeholder.R
  tests/
    testthat/
      test_seurat_input.R
      test_pooling_na.R
      test_empty_pool.R
      test_pseudobulk_order.R
      test_safe_scale.R
      test_boltzmann_and.R
      test_q95_continuous.R
      test_confidence_nonnegative.R
  inst/
    extdata/
      toy_gpr.tsv
      toy_pool_meta.tsv
  scripts/
    01_run_layer1_from_seurat_v4.R
    02_diagnostics_report.R
  output/
    pool/
    layer1/
    confidence/
    diagnostics/
```

---

## 2. Seurat v4 输入规格

### 2.1 必需输入

Seurat v4 对象必须包含：

```text
RNA assay: usually "RNA" or "SCT"
ATAC assay: usually "ATAC" / "peaks"
metadata:
  sample_id
  cell_type or major_cell_type
optional metadata:
  condition
  batch
  local_state / subcluster / epithelial_substate
  ATAC QC: TSS.enrichment, FRiP, nucleosome_signal, peak_region_fragments
```

输入对象已经注释好，因此 RegCompassR 不负责重新注释 cell type。

### 2.2 推荐输入

```text
gpr_table.tsv:
  reaction_id
  gpr
  pathway
  compartment
  is_transport
  is_exchange

peak_gene_links.tsv:
  peak_id
  gene
  weight
  link_source
  distance_to_TSS

motif_deviation matrix:
  TF/motif × cell or TF/motif × pool
```

第一版可以只要求 `Seurat object + gpr_table.tsv`，ATAC confidence 可作为 v0.4 模块加入。

---

## 3. 开发路线图

## v0.1：Seurat v4 输入与 metadata 验证

目标：

```text
确认 Seurat v4 对象结构正确
抽取 RNA counts / ATAC counts / metadata
检查 sample_id、cell_type、condition、state_col
```

核心函数：

```r
rc_validate_seurat_v4()
rc_extract_seurat_v4()
rc_check_metadata()
rc_write_input_summary()
```

### 3.1 输入检查规则

必须检查：

```text
1. object 是否继承 Seurat
2. RNA assay 是否存在
3. ATAC assay 是否存在
4. RNA 和 ATAC 是否有相同 cell barcode
5. sample_col 是否存在
6. celltype_col 是否存在
7. condition_col 若指定则必须存在
8. state_col 若指定则必须存在
9. sample × cell_type 的细胞数分布
10. condition 与 batch 是否明显混杂
```

### 3.2 Seurat v4 数据抽取

```r
rc_get_assay_counts <- function(object, assay) {
  SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
}
```

说明：

- 针对 Seurat v4，`slot="counts"` 是主要接口。
- 不从 `scale.data` 或 imputed matrix 计算 GPR capacity。
- 若将来兼容 Seurat v5，再增加 `layer="counts"` 分支。

---

## v0.2：Sample-aware micropooling

## 4. Pool 定义

主定义：

```text
pool = sample_id × condition × cell_type × local_state
```

若无处理条件：

```text
pool = sample_id × cell_type × local_state
```

若无 `local_state`：

```text
pool = sample_id × condition × cell_type
```

### 4.1 关于 local_state 的修正原则

Biomni 批判指出 `local_state` 生物学地位不清，这一点成立。因此：

1. `state_col` 必须由用户显式指定，不能自动使用任意 `seurat_clusters`。
2. 如果 `state_col = seurat_clusters`，必须在 diagnostics 中记录 resolution / clustering source。
3. 如果 `state_col` 来自 WNN cluster，应标记 `state_source = WNN`，提示后续 ATAC confidence 可能与 state 定义不完全独立。
4. 第一版默认允许 `state_col = NULL`，只按 `sample × cell_type` 建 pool。
5. 若要研究上皮内部状态，推荐使用人工注释的 `epithelial_substate`，而不是任意分辨率聚类。

### 4.2 Pool 大小参数

默认：

```text
target_pool_size = 80
min_pool_size = 30
min_group_size = 30
```

但必须声明：

```text
这些不是生物学定律，只是计算降噪与状态分辨率的折中。
```

必须输出 sensitivity：

```text
target_pool_size = 50, 80, 120
min_pool_size = 20, 30
```

### 4.3 对单一 pool group 的处理

如果某个 `sample × cell_type × state` group 只有一个 pool：

```text
保留该 pool
但标记 no_within_group_pool_replicate = TRUE
不得用它估计 within-sample pool variation
```

注意：pool 本来不是生物学重复。一个 group 只有一个 pool 不等于不能分析，但不能把它用于 pool-level 方差估计。

### 4.4 随机性控制

`rc_make_pools()` 使用随机打乱 cell，因此必须：

```text
seed 必须记录
pool_seed_replicates 可选
Layer 1 sensitivity 可运行多个 seed
```

推荐：

```text
pool_seed_replicates = 1 for default
pool_seed_replicates = 5 for robustness check
```

---

## 5. Pooling 实现细节

### 5.1 NA 分组处理

严禁 NA 进入 grouping key。

```r
rc_drop_na_grouping <- function(meta, group_cols) {
  bad <- rowSums(is.na(meta[, group_cols, drop = FALSE])) > 0
  if (any(bad)) {
    warning(sum(bad), " cells removed due to NA in grouping columns")
  }
  meta[!bad, , drop = FALSE]
}
```

### 5.2 Pool 构建伪代码

```r
rc_make_pools <- function(
  meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = NULL,
  state_col = NULL,
  target_size = 80,
  min_pool_size = 30,
  min_group_size = 30,
  seed = 1
) {
  set.seed(seed)
  group_cols <- c(sample_col, condition_col, celltype_col, state_col)
  group_cols <- group_cols[!is.null(group_cols)]

  meta <- rc_drop_na_grouping(meta, group_cols)
  meta$cell_id <- rownames(meta)

  keys <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, sep = "|")
  groups <- split(meta$cell_id, keys)

  out <- list()
  k <- 1L

  for (nm in names(groups)) {
    cells <- groups[[nm]]
    n <- length(cells)

    if (n < min_group_size) {
      out[[k]] <- data.frame(
        group_id = nm,
        pool_id = NA_character_,
        cell_id = cells,
        skipped = TRUE,
        skip_reason = "group_below_min_group_size"
      )
      k <- k + 1L
      next
    }

    n_pool <- max(1L, floor(n / target_size))
    cells <- sample(cells)
    pool_assign <- rep(seq_len(n_pool), length.out = n)

    for (j in seq_len(n_pool)) {
      cc <- cells[pool_assign == j]
      out[[k]] <- data.frame(
        group_id = nm,
        pool_id = paste0("pool_", k),
        cell_id = cc,
        skipped = FALSE,
        low_power_pool = length(cc) < min_pool_size,
        no_within_group_pool_replicate = n_pool == 1L
      )
      k <- k + 1L
    }
  }

  do.call(rbind, out)
}
```

---

## v0.3：Raw-count pseudobulk 与 pool-level normalization

## 6. Pseudobulk 顺序

必须使用：

```text
raw counts → pool sum → pool-level normalization
```

不得使用：

```text
mean(cell-level Pearson residual)
mean(cell-level TF-IDF)
mean(imputed expression)
```

原因：

```text
mean(PearsonResidual) ≠ PearsonResidual(pseudobulk counts)
mean(TFIDF) ≠ TFIDF(pseudobulk counts)
```

### 6.1 RNA pseudobulk

```r
rc_pseudobulk_counts <- function(counts, pool_map, fun = "sum") {
  pool_map <- pool_map[!is.na(pool_map$pool_id) & !pool_map$skipped, ]
  pool_ids <- unique(pool_map$pool_id)

  res <- lapply(pool_ids, function(pid) {
    cells <- pool_map$cell_id[pool_map$pool_id == pid]
    x <- counts[, cells, drop = FALSE]
    if (fun == "sum") Matrix::rowSums(x) else Matrix::rowMeans(x)
  })

  out <- do.call(cbind, res)
  colnames(out) <- pool_ids
  out
}
```

### 6.2 空 pool 处理

必须在 normalization 前过滤 library size 为 0 的 pool。

```r
rc_filter_empty_pools <- function(pb_counts, pool_meta) {
  lib <- Matrix::colSums(pb_counts)
  keep <- lib > 0
  if (any(!keep)) {
    warning(sum(!keep), " empty pools removed before normalization")
  }
  list(
    counts = pb_counts[, keep, drop = FALSE],
    pool_meta = pool_meta[match(colnames(pb_counts)[keep], pool_meta$pool_id), , drop = FALSE]
  )
}
```

### 6.3 RNA pool-level normalization

第一版使用 logCPM：

\[
X^{RNA}_{g,p} = \log_2\left(1 + \frac{count_{g,p}}{\sum_g count_{g,p}} \times 10^6\right)
\]

```r
rc_logcpm <- function(pb_counts, scale_factor = 1e6) {
  lib <- Matrix::colSums(pb_counts)
  if (any(lib <= 0)) stop("Empty pools detected. Run rc_filter_empty_pools() first.")
  norm <- t(t(pb_counts) / lib) * scale_factor
  log1p(norm) / log(2)
}
```

说明：

- 该值是 `log2(CPM + 1)`。
- 对低表达基因有压缩作用；这需要在 diagnostics 中报告 low-expression gene fraction。
- 后续可增加 `scran` size factor 或 DESeq2-style normalization 作为 sensitivity。

---

## 7. ATAC pool-level normalization：不以 TF-IDF 作为代谢 confidence 主输入

Biomni 批判指出 pool-level TF-IDF 会惩罚所有 pool 中广泛开放的 peaks，这对 housekeeping/metabolic gene promoter 可能不合适。该问题成立。

因此调整如下：

### 7.1 ATAC regulatory support 主输入

第一版默认使用 pool-level accessibility logCPM：

\[
X^{ATAC}_{e,p} = \log_2\left(1 + \frac{count_{e,p}}{\sum_e count_{e,p}} \times 10^6\right)
\]

然后在 stratum 内转 percentile：

\[
P^{ATAC}_{e,p} = percentile_{p \in stratum}(X^{ATAC}_{e,p})
\]

### 7.2 TF-IDF 的位置

TF-IDF 只用于：

```text
ATAC LSI / embedding / optional peak specificity diagnostics
```

不作为 `LinkConf` 或代谢基因 regulatory support 的默认输入。

### 7.3 全零 peak 过滤

在 pool-level ATAC normalization 前过滤：

```text
peak detected in at least min_pools peaks
default min_pools = 3
```

---

## v0.4：Gene score 与 GPR Layer 1 capacity

## 8. Gene score

### 8.1 Safe robust scale

Biomni 批判指出 `MAD(constant=1)` 与 `IQR/1.35` 量纲不一致，该问题成立。因此统一使用 σ 估计：

\[
scale_g = \max(MAD_g^{\sigma}, IQR_g/1.349, min\_scale)
\]

其中：

\[
MAD_g^{\sigma} = median(|x - median(x)|) \times 1.4826
\]

默认：

```text
min_scale = 0.05
z_clip = 6
```

```r
rc_safe_scale <- function(x, min_scale = 0.05) {
  mad_sigma <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  iqr_sigma <- stats::IQR(x, na.rm = TRUE) / 1.349
  max(mad_sigma, iqr_sigma, min_scale, na.rm = TRUE)
}

rc_gene_zscore <- function(X, min_scale = 0.05, z_clip = 6) {
  z <- X
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ])
    med <- stats::median(x, na.rm = TRUE)
    sc <- rc_safe_scale(x, min_scale = min_scale)
    zi <- (x - med) / sc
    zi <- pmax(pmin(zi, z_clip), -z_clip)
    z[i, ] <- zi
  }
  z
}

rc_sigmoid <- function(z) 1 / (1 + exp(-z))
```

### 8.2 Detection rate 不直接修正 gene score

Biomni 批判指出：

```text
q_gp 是 cell-level detection rate
s_gp 来自 pool-level logCPM
二者直接线性混合存在尺度不匹配
```

因此主分析删除：

\[
\tilde{s}_{g,p} = q_{g,p}s'_{g,p} + (1-q_{g,p})m_{g,c}
\]

主分析改为：

\[
s_{g,p} = \sigma(z_{g,p})
\]

detection rate 只用于：

```text
gene_confidence
low_detection_flag
dropout_sensitive_flag
```

可选 sensitivity：

```text
只在 q_high pools 中估计 baseline m_gc
但不作为默认
```

### 8.3 Promiscuity correction

保留三种模式：

\[
s'_{g,p} =
\begin{cases}
s_{g,p}, & none \\
s_{g,p}/\sqrt{N_{rxn}(g)}, & sqrt \\
s_{g,p}/N_{rxn}(g), & linear
\end{cases}
\]

默认：

```text
promiscuity = "sqrt"
```

但必须报告：

```text
none / sqrt / linear sensitivity
```

---

## 9. GPR parser

第一版支持标准 GPR：

```text
(g1 and g2) or g3
g1 or g2
g1 and g2
g1
```

复杂嵌套建议使用正式 parser，不要长期依赖简单 `strsplit`。第一版可先接受预解析表：

```text
reaction_id
and_group_id
gene
```

推荐内部格式：

| reaction_id | and_group_id | gene |
|---|---:|---|
| R1 | 1 | GENE1 |
| R1 | 1 | GENE2 |
| R1 | 2 | GENE3 |

这样可避免字符串解析错误。

---

## 10. AND aggregation：Boltzmann-weighted average 与 τ

### 10.1 算子定义

\[
w_{g,p}=\frac{\exp(-s'_{g,p}/\tau)}{\sum_{h\in G_{r,k}}\exp(-s'_{h,p}/\tau)}
\]

\[
C_k(r,p)=\sum_{g\in G_{r,k}}w_{g,p}s'_{g,p}
\]

性质：

\[
\min(s')\le C_k\le mean(s')
\]

### 10.2 τ=0.08 的重新定位

Biomni 批判指出：在 sigmoid score \([0,1]\) 中，τ=0.08 常接近 hard min。该问题成立。

因此：

```text
τ=0.08 不再作为“soft AND”默认解释；
τ=0.08 = strict bottleneck mode；
τ=0.20 或 0.30 = softer AND mode；
hard min = reference baseline；
mean = upper sensitivity。
```

### 10.3 推荐默认

v0.4 默认运行：

```text
AND methods:
  min
  boltzmann_tau_0.08
  boltzmann_tau_0.20
  mean
```

主报告可选择：

```text
default_and = boltzmann_tau_0.20
```

但必须输出：

```text
and_method_sensitivity
tau_sensitive_flag
```

### 10.4 缺失亚基处理

如果 AND group 中某些 gene 缺失：

```text
missing_subunit_fraction = missing / total
if missing_subunit_fraction > 0:
  capacity still computed using available genes
  confidence is downgraded
  missing_subunit_flag = TRUE
```

严禁静默丢弃 NA 而不报告。

---

## 11. OR aggregation

\[
C_{raw}(r,p)=\sum_{k=1}^{K_r} C_k(r,p)
\]

解释：

```text
OR 表示 isoenzyme / alternative enzymatic route
加和保留 isoenzyme 累积容量
```

输出：

```text
C_raw
n_isoenzyme_groups
has_isoenzyme
has_multisubunit
missing_gene_fraction
```

---

## v0.5：Q95 continuous calibration

## 12. 连续 Q95 shrinkage

Biomni 批判指出 `n<20` 硬切换到 global Q 会破坏连续性，这一点成立。因此修正为：

\[
Q_r = \rho_n Q_{r,stratum} + (1-\rho_n)Q_{r,global}
\]

\[
\rho_n=\frac{n}{n+n_0}
\]

默认：

```text
n0 = 80
q = 0.95
```

适用于所有 \(n\ge1\)。

### 12.1 低功效标记

不切换公式，只打标记：

```text
n < 5: q95_very_low_power
n < 20: q95_low_power
n < 100: q95_moderate_power
n ≥ 400: q95_high_power
```

### 12.2 Q95 bootstrap

默认：

```text
B = 500 for regular mode
B = 1000 for accurate mode
```

输出：

```text
q_stratum
q_global
rho_n
q_shrink
q95_ci_low
q95_ci_high
q95_ci_width
q95_unstable_flag
```

### 12.3 C_rel

\[
C_{rel}(r,p)=\min\left(1,\frac{C_{raw}(r,p)}{Q_r+\epsilon}\right)
\]

注意：

```text
C_rel 是相对 capacity potential
不是绝对酶活性
不是 flux bound
```

---

## v0.6：Multiome confidence 重构

## 13. Stratum-wise percentile

所有 percentile 默认在同一 stratum 内计算：

```text
stratum = cell_type
optional: cell_type × condition
optional: cell_type × local_state
```

不建议跨全部 cell type 计算 percentile，否则细胞类型身份会主导 TFAct / ATAC support。

---

## 14. RNA-ATAC concordance

### 14.1 原始 concordance

\[
Concord_{RA}(g,p)=1-|P_{RNA}(g,p)-P_{ATAC}(g,p)|
\]

### 14.2 离散 null baseline

对 n 个 pool，独立随机 rank 的期望 concordance 为：

\[
E_{null}(n)=\frac{2}{3}+\frac{1}{3n^2}
\]

因此：

\[
Concord^{norm}_{RA}(g,p)=
\max\left(0,\frac{Concord_{RA}(g,p)-E_{null}(n)}
{1-E_{null}(n)}\right)
\]

当 n 很小时，必须使用该离散 correction，而不是固定 \(2/3\)。

---

## 15. Fisher shrinkage

Spearman 相关：

\[
\rho_{RA}(g)=Spearman(X^{RNA}_{g,\cdot}, X^{ATAC}_{g,\cdot})
\]

clip：

\[
\rho^{clip}=\min(0.999,\max(-0.999,\rho))
\]

Fisher z：

\[
z=\operatorname{arctanh}(\rho^{clip})
\]

收缩：

\[
\lambda=\max\left(0,\frac{n-3}{n-3+n_0}\right)
\]

\[
\rho^{shrink}=\tanh(\lambda z)
\]

默认：

```text
n0 = 30
n < 4: do not compute correlation, rho_shrink = 0
n < 10: low_correlation_power_flag = TRUE
```

说明：

```text
Fisher z 对 Spearman 是近似方法；
小 n 时结果仅作为 diagnostics，不作为强证据。
```

---

## 16. Gene confidence：非负重构

Biomni 批判指出 gene confidence 可能为负，这一点成立。修正如下。

### 16.1 正相关支持

\[
Rel^+_{RA}(g)=\max(0,\rho^{shrink}_{RA}(g))
\]

负相关不作为支持，而是单独记录：

\[
Discord_{RA}(g)=|\min(0,\rho^{shrink}_{RA}(g))|
\]

### 16.2 Confidence 公式

\[
Conf(g,p)=
0.25\cdot Concord^{norm}_{RA}(g,p)Rel^+_{RA}(g)
+
0.15\cdot Concord^{norm}_{RT}(g,p)Rel^+_{RT}(g)
+
0.20\cdot Det_{RNA}(g,p)
+
0.15\cdot LinkConf(g,p)
+
0.15\cdot QC(p)
+
0.10\cdot GPRGeneObserved(g,p)
\]

所有分量均在 \([0,1]\)，因此：

\[
Conf(g,p)\in[0,1]
\]

### 16.3 TFAct

默认：

\[
TFAct(f,p)=\frac{1}{2}\left(P_{RNA}(f,p)+P_{motif}(f,p)\right)
\]

其中：

```text
P_RNA = TF RNA expression percentile within stratum
P_motif = chromVAR deviation percentile within stratum
```

说明：

```text
TFAct 是 proxy，不是真实 TF activity。
对 SREBP、HIF1A、MYC、ChREBP 等代谢 TF，mRNA/motif proxy 可能不足。
```

### 16.4 LinkConf

使用 pool-level ATAC normalized accessibility，而非 TF-IDF：

\[
LinkConf(g,p)=\sum_{e\in Enh(g)} w'_{e,g}P^{ATAC}_{e,p}
\]

权重归一化：

\[
w'_{e,g}=\frac{\max(0,w_{e,g})}{\sum_{e}\max(0,w_{e,g})}
\]

若 peak-gene weight 有负值：

```text
negative links are stored as repressive_link_flag
not used as positive support
```

---

## 17. Reaction confidence

\[
RConf(r,p)=weighted\_median\{Conf(g,p): g\in GPR(r)\}
\]

默认权重：

```text
AND group 内：低 capacity gene 权重更高
OR group 间：group capacity 权重
```

第一版可简化为：

```text
RConf(r,p) = median Conf(g,p) over GPR genes
```

但必须报告：

```text
missing_gpr_gene_fraction
low_confidence_reaction_flag
```

---

## v0.7：Layer 1 输出与 sample-level 汇总

## 18. 第一版主输出

第一版不输出 flux，只输出：

```text
reaction_capacity_raw
reaction_capacity_relative
reaction_confidence
q95_diagnostics
multiome_confidence
pool_metadata
sample_level_summary
```

## 19. Sample-level aggregation

\[
Y_{s,c,r}=median_{p\in(s,c)} C_{rel}(r,p)
\]

同时输出：

```text
n_pools_used
n_cells_used
single_pool_group_flag
low_power_group_flag
```

如果有 condition：

\[
Y_{s,t,c,r}=median_{p\in(s,t,c)} C_{rel}(r,p)
\]

第一版不建议自动做复杂统计模型，只输出可供用户建模的 sample-level matrix。

---

## v0.8：并行与 checkpoint

## 20. 并行原则

Linux R 使用 `BiocParallel`：

```r
BPPARAM = BiocParallel::MulticoreParam(workers = 8)
```

原则：

```text
1. 不把整个 Seurat object 传入 worker
2. 只传递必要的 matrix block 或 reaction chunk
3. 不嵌套并行
4. 每个耗时步骤支持 checkpoint
5. 支持 resume
```

主要并行步骤：

```text
GPR reaction capacity by reaction chunks
Q95 bootstrap by reaction
seed sensitivity by seed
future QP by pool
```

---

## v0.9：未来 QP 层的重新定位

Biomni 批判指出 capacity 作为 hard flux bound 假设过强，该问题成立。

因此 QP 层只作为后续模块，且重新定位为：

```text
selected network-constrained feasibility
```

不是：

```text
true flux estimation
```

默认 QP 设计：

```text
C_rel → objective soft penalty
not default hard upper bound
```

hard capacity bound 仅用于 sensitivity：

```text
Mode B: hard-bound sensitivity
```

Transport/exchange reaction：

```text
不因无 GPR 且 C_rel=1 而自动解释为高 confidence
只有用户指定、GEM 注释为关键 transport/exchange 或与 selected pathway 有关时进入 selected QP
```

---

## 21. 版本化开发计划

## v0.1：Seurat v4 输入

实现：

```text
rc_validate_seurat_v4()
rc_extract_seurat_v4()
rc_check_metadata()
rc_write_input_summary()
```

完成标准：

```text
读取已注释 Seurat v4 multiome
检查 RNA/ATAC cells 一致
检查 sample_id/cell_type/condition/state_col
```

---

## v0.2：Pooling + pseudobulk

实现：

```text
rc_make_pools()
rc_pseudobulk_counts()
rc_filter_empty_pools()
rc_logcpm()
rc_pool_metadata()
```

完成标准：

```text
无跨 sample pool
无 NA group
无空 pool normalization
可输出 pool_metadata.tsv
```

---

## v0.3：Layer 1 GPR capacity

实现：

```text
rc_parse_gpr_table()
rc_gene_zscore()
rc_gene_score()
rc_promiscuity()
rc_and_capacity()
rc_or_capacity()
rc_q95_shrink()
rc_layer1_capacity()
```

完成标准：

```text
toy GPR 测试通过
AND sensitivity 输出
promiscuity sensitivity 输出
C_raw/C_rel 输出
```

---

## v0.4：Multiome confidence

实现：

```text
rc_atac_pool_logcpm()
rc_percentile_by_stratum()
rc_concordance_null_correct()
rc_fisher_shrink()
rc_link_confidence()
rc_gene_confidence()
rc_reaction_confidence()
```

完成标准：

```text
confidence 非负且 ≤1
负相关单独输出 discordance
TFAct 用 stratum-wise percentile
ATAC 不默认使用 TF-IDF
```

---

## v0.5：Diagnostics report

实现：

```text
rc_diagnostics_pool()
rc_diagnostics_gpr()
rc_diagnostics_q95()
rc_diagnostics_confidence()
rc_write_report_md()
```

必须报告：

```text
pool size distribution
single_pool_group fraction
q95 low-power fraction
tau sensitivity
promiscuity sensitivity
confidence distribution
discordance fraction
missing GPR genes
```

---

## v0.6：Sample-level summary

实现：

```text
rc_sample_aggregate()
rc_export_sample_matrix()
rc_export_long_table()
```

输出：

```text
sample × cell_type × reaction matrix
pool-level long table
reaction metadata
diagnostics table
```

---

## v0.7+：Selected QP prototype

仅在 Layer 1 与 confidence 稳定后开发。

实现前置条件：

```text
1. Layer 1 seed sensitivity acceptable
2. Q95 unstable fraction low
3. GPR missing fraction reported
4. confidence not dominated by QC/detection
5. sample-level results biologically interpretable
```

---

## 22. 最小可运行流程

```r
library(Seurat)
library(Signac)
library(Matrix)
library(data.table)
library(BiocParallel)

seu <- readRDS("annotated_multiome_seurat_v4.rds")

rc_validate_seurat_v4(
  object = seu,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  state_col = NULL
)

inp <- rc_extract_seurat_v4(
  object = seu,
  rna_assay = "RNA",
  atac_assay = "ATAC"
)

pool_map <- rc_make_pools(
  meta = inp$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  state_col = NULL,
  target_size = 80,
  min_pool_size = 30,
  min_group_size = 30,
  seed = 1
)

rna_pb <- rc_pseudobulk_counts(inp$rna_counts, pool_map)
filtered <- rc_filter_empty_pools(rna_pb, rc_build_pool_metadata(pool_map, inp$meta))
rna_pb <- filtered$counts
pool_meta <- filtered$pool_meta

rna_logcpm <- rc_logcpm(rna_pb)

gpr <- data.table::fread("gpr_table_long.tsv")

layer1 <- rc_layer1_capacity(
  rna_pool = rna_logcpm,
  gpr_table = gpr,
  pool_meta = pool_meta,
  and_methods = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"),
  promiscuity_modes = c("none", "sqrt", "linear"),
  q95_n0 = 80,
  BPPARAM = BiocParallel::MulticoreParam(workers = 8)
)

data.table::fwrite(layer1$capacity_long, "output/layer1_capacity.tsv")
data.table::fwrite(layer1$q95_diagnostics, "output/q95_diagnostics.tsv")
data.table::fwrite(layer1$gpr_diagnostics, "output/gpr_diagnostics.tsv")
```

---

## 23. 必须加入的单元测试

| 测试 | 目的 |
|---|---|
| `test_pooling_no_cross_sample` | 确认 pool 不跨 sample |
| `test_pooling_na_excluded` | 确认 NA grouping cells 被排除 |
| `test_empty_pool_removed` | 防止空 pool 被 logCPM 救援 |
| `test_pseudobulk_raw_before_norm` | 确认 raw-count pseudobulk 顺序 |
| `test_safe_scale_sigma` | 确认 MAD/IQR 量纲一致 |
| `test_boltzmann_range` | 确认 min ≤ AND ≤ mean |
| `test_tau_008_min_like` | 确认 τ=0.08 接近 min 并正确标记 |
| `test_q95_no_hard_switch` | 确认所有 n 连续 shrinkage |
| `test_confidence_nonnegative` | 确认 Conf ∈ [0,1] |
| `test_negative_correlation_discordance` | 确认负相关进入 discordance 不降低为负权重 |
| `test_concordance_null_discrete` | 确认小 n 使用离散 null baseline |
| `test_detection_not_modify_score` | 确认 detection rate 不直接修正 gene score |

---

## 24. 最终结论

根据 Biomni 批判报告，本方案的核心修订是：

```text
1. 第一版只做 Layer 1 capacity，不急于做 QP/FVA。
2. 输入明确限定为已注释 Seurat v4 multiome 对象。
3. Pseudobulk 必须 raw counts 先聚合，再 pool-level normalization。
4. ATAC confidence 默认不使用 TF-IDF，而使用 pool-level normalized accessibility。
5. Detection rate 不再直接修正 gene score，只作为 confidence/diagnostics。
6. AND 聚合不再固定解释 τ=0.08 为 soft；τ=0.08 作为 strict-min-like sensitivity。
7. Q95 对所有 n 使用连续 shrinkage，不做硬切换。
8. Gene confidence 必须非负；负相关作为 discordance flag。
9. QP 只作为后续 selected feasibility 模块，不作为 true flux 或主版本目标。
```

该路线比直接实现完整 RegCompass-Multiome 更稳健。当前最合理的 MVP 是：

```text
Seurat v4 annotated multiome
→ sample-aware raw-count pseudobulk
→ pool-level RNA/ATAC normalization
→ GPR Layer 1 reaction capacity
→ nonnegative multiome confidence
→ diagnostics + sample-level summary
```

只有当 Layer 1 的 seed sensitivity、GPR coverage、Q95 stability 和 confidence diagnostics 均通过后，才应进入 selected Human-GEM QP 阶段。
