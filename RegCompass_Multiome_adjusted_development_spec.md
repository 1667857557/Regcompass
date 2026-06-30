# RegCompass-Multiome：单细胞同胞 RNA+ATAC 多组学代谢分析工具的严谨审查与开发规格

> 输入：Biomni 输出的 RegCompass-Multiome 方案讨论与后续参数/数学验证报告。  
> 目标：判断报告是否正确，并在保留正确计算验证细节的基础上，形成一份可用于工具开发的完整 Markdown 规格。  
> 最终定位：RegCompass-Multiome 不是真实 flux 估计器，而是 **multiome-supported, GPR-aware, sample-aware 的 reaction capacity potential 与 selected network-constrained feasibility 分析框架**。

---

## 0. 总体判断

Biomni 讨论报告的主体方向是正确的，尤其是以下判断：

1. RNA 和 ATAC 不能直接解释为真实酶活性或真实代谢通量。
2. ATAC/motif 应作为 regulatory support / confidence，而不是 flux bound。
3. sample-aware micropooling 是合理的降噪和计算加速策略，但 pool 不是独立生物学重复。
4. 对 GPR 中 AND/OR 关系进行显式建模是必要的。
5. 全量 reaction × pool 的 exact QP 求解不可行，应使用 Layer 1 全量快速评分 + Layer 3 selected exact feasibility。
6. Regulator ranking 只能称为 candidate prioritization，不能称 causal driver discovery。

但报告仍需修正或补充：

| 问题 | 原报告状态 | 修正结论 |
|---|---|---|
| AND “softmin” 命名 | 称为 softmin | 应称为 **Boltzmann-weighted average / soft-minimum via weighted average**，不是 LogSumExp softmin |
| METAFlux 对比 | 描述较浅 | 必须明确 METAFlux AND = hard min，OR = sum，promiscuity = 1/N，RegCompass 与其不同 |
| QP 使用 | 已识别双重惩罚风险 | 主分析采用 **soft-expression penalty**，hard capacity bound 仅作 sensitivity |
| Q95 校准 | 已提出 n 阈值 | 必须增加 bootstrap CI；absolute floor 只能作 flag，不能作硬过滤 |
| dropout correction | 已识别全局 baseline 风险 | 必须改为 cell-type/substate-specific baseline |
| FVA | 已识别 biomass objective 约束错误 | ATP mode 用 objective tolerance；biomass mode 用 biomass optimum 保持约束 |
| 可逆性/区室/热力学 | 原报告不足 | 必须加入 reversibility、futile cycle、compartment、flux consistency、loop diagnostics |
| 计算成本 | 原报告不足 | 必须在工具中预估 QP 数量、运行时间、并行计划与 checkpoint |

---

## 1. 参考工具算法定位与对 RegCompass 的约束

### 1.1 scMetabolism

**算法定位**：基于 gene-set scoring 的单细胞代谢通路活性评分工具，支持 VISION、AUCell、ssGSEA、GSVA；数据库包括 KEGG 和 Reactome metabolic pathway。  
**不包含**：GEM、GPR、stoichiometric constraints、FBA/QP、reaction-level feasibility。

**对 RegCompass 的启示**：

- 保留一个快速 Layer 0：AUCell/UCell/scMetabolism-like pathway prescreening。
- Layer 0 只用于可视化和候选 pathway/reaction 筛选，不能作为 reaction-level 或 flux 证据。

---

### 1.2 scFEA

**算法定位**：将代谢网络压缩为 factor graph / module network，用 GNN 从 scRNA-seq 预测 cell-wise module flux-like profile。其优势是原生 cell-wise 输出和非线性表达映射，但其 module-level 结果低于 genome-scale reaction 分辨率。

**对 RegCompass 的启示**：

- 可借鉴 network consistency / balance loss 的思想。
- 不应照搬 GNN 黑箱；RegCompass 的核心优势应是 GPR-explicit、reaction-level、diagnostic-rich。
- 不应声称优于 scFEA 的 cell-wise 分辨率；RegCompass 是 pool-level main analysis + cell-wise projection。

---

### 1.3 COMPASS

**算法定位**：使用 scRNA-seq 与 FBA/LP 框架计算 reaction penalty / consistency。COMPASS 支持 micropooling、selected reactions、penalty diffusion；其 `--and-function` 支持 `min/median/mean`，默认实现需随版本和命令记录。COMPASS 文档明确指出，使用 micropools/metacells 后通常不需要 information sharing，`lambda` 应设为 0。

**对 RegCompass 的启示**：

- 继承 reaction-level 解释性和 selected exact scoring 策略。
- micropooling 后默认关闭 KNN diffusion。
- 对照 COMPASS 时，不应直接比较数值尺度；只能比较 rank concordance 或 direction consistency。

---

### 1.4 METAFlux

**算法定位**：以 Human1/Human-GEM 为底模，使用 GPR-derived MRAS、nutrient profile 和 QP-FBA 预测 reaction flux scores。METAFlux 公开文献描述其可输出 13,082 reaction flux scores，且强调 nutrient-aware modeling。

METAFlux MRAS 的关键逻辑可抽象为：

- OR / isoenzyme：

\[
MRAS_{r} = \sum_{k=1}^{K_r} \frac{Enzyme_k}{w_k}
\]

- AND / complex：

\[
MRAS_{r} = \min_k \left(\frac{Enzyme_k}{w_k}\right)
\]

- promiscuity correction：

\[
w_k = N_{rxn}(gene_k)
\]

即 METAFlux 对 promiscuous gene 采用近似 \(1/N_{rxn}\) 的线性惩罚。

**RegCompass 与 METAFlux 的重要差异**：

| 维度 | METAFlux | RegCompass-Multiome |
|---|---|---|
| 输入 | bulk/scRNA | paired scRNA+scATAC multiome |
| AND | hard min | Boltzmann-weighted average，抗 dropout |
| OR | additive sum | additive raw sum + Q95 calibration |
| promiscuity | \(1/N_{rxn}\) | 默认 \(1/\sqrt{N_{rxn}}\)，并做 sensitivity |
| expression→QP | 更接近 hard bound / MRAS constraint | 默认 soft penalty，hard bound 仅 sensitivity |
| objective | biomass/nutrient-aware flux prediction | ATP maintenance baseline + selected demand feasibility |
| TME exchange | 可建模 community exchange | 默认删除，除非有空间/代谢物交换证据 |

---

## 2. 生物学边界

RegCompass-Multiome 的推断链为：

\[
ATAC/motif \rightarrow RNA \rightarrow protein \rightarrow active\ enzyme \rightarrow catalytic\ capacity \rightarrow flux
\]

实际数据只覆盖前两层。因此：

| 输出术语 | 可以使用 | 不应使用 |
|---|---|---|
| gene score | gene expression-derived capacity score | enzyme activity |
| reaction score | reaction capacity potential | reaction flux |
| QP result | network-constrained feasibility | true intracellular flux |
| ATAC/motif | regulatory support confidence | reaction activation proof |
| TF ranking | candidate regulator / mediator | causal driver |
| cell score | projected visualization score | independently solved single-cell flux |

### 2.1 ATAC 时间尺度限制

ATAC peak、motif accessibility、enhancer-gene link 反映的是染色质调控潜势，常处于小时到天尺度；flux 变化可在秒到分钟尺度发生。因此 ATAC 不能进入 flux bound，只能作为：

```text
confidence modifier
regulatory support annotation
candidate regulator ranking evidence
state stratification evidence
```

### 2.2 乳腺癌/肿瘤数据特殊注意

对于恶性上皮细胞或乳腺癌细胞：

- 应控制 CNV dosage / malignant clone / subtype / cell cycle。
- EGFR/HER2/ERBB pathway 分析应使用 pathway activity，而非单基因 EGFR 表达。
- metabolic rewiring 的解释必须分清 basal-like、EMT-like、cycling、stress-like、HER2/EGFR-high 等状态混杂。

---

## 3. 最终工具目标

RegCompass-Multiome 应输出：

1. gene-level RNA/ATAC/motif support
2. GPR-aware reaction capacity potential
3. selected network-constrained reaction feasibility
4. pathway / meta-reaction potential
5. regulator candidate ranking
6. sample-aware differential results
7. cell-wise projected visualization score
8. diagnostics report

明确不输出：

```text
true flux
true enzyme activity
causal regulatory proof
single-cell independently solved flux
TME exchange flux without direct evidence
thermodynamically validated flux unless explicit thermodynamic constraints are used
```

---

## 4. 输入与输出规格

### 4.1 必需输入

```text
RNA count matrix: genes × cells
ATAC peak matrix: peaks × cells
cell metadata: sample, condition, batch, cell_type, QC
joint embedding: WNN / MOFA / totalVI-like representation
peak-gene links: enhancer → gene with weights and source
motif deviation / TF activity matrix
Human-GEM / Human1 model with GPR rules
sample-level biological replicate labels
```

### 4.2 推荐输入

```text
nutrient / medium profile
malignant CNV dosage
cell cycle score
protein / Ribo-seq / phosphoproteomics
matched metabolomics
13C isotope tracing
Seahorse OCR/ECAR
perturbation multiome
```

### 4.3 输出文件

```text
01_pool_metadata.tsv
02_gene_multiome_support.tsv
03_reaction_capacity_L1.tsv
04_qp_baseline_diagnostics.tsv
05_reaction_feasibility_L3.tsv
06_limited_FVA.tsv
07_pathway_meta_reaction_scores.tsv
08_sample_level_differential_results.tsv
09_regulator_candidate_ranking.tsv
10_cell_projected_scores.tsv
11_validation_sensitivity_summary.tsv
12_diagnostics_report.md
```

---

## 5. 计算模块

## Module 0：QC、normalization 与 joint embedding

### RNA QC

必须记录：

```text
nUMI
nGene
mitochondrial fraction
ribosomal fraction
ambient contamination estimate
doublet status
```

### ATAC QC

必须记录：

```text
TSS enrichment
FRiP
nucleosome signal
peak count
blacklist ratio
```

### Normalization

RNA 主体分析建议使用 sample-aware Pearson residual 或 SCTransform residual：

\[
Z_{RNA}(g,c)=PearsonResidual(g,c)
\]

ATAC 使用 TF-IDF + LSI：

\[
A(e,c)=TFIDF(e,c)
\]

joint embedding 可使用 WNN：

```text
RNA PCA + ATAC LSI → WNN graph → local state clustering
```

**实现注意**：

- corrected embedding 可用于聚类和 local state。
- detection rate、UMI、QC 必须从 raw / normalized-but-not-imputed 数据计算。
- imputation 不得用于 GPR capacity 主公式。

---

## Module 1：Sample-aware micropooling

### Pool 定义

\[
pool \subset sample \times condition \times cell\_type \times local\_state
\]

若没有 treatment contrast，可省略 condition：

\[
pool \subset sample \times cell\_type \times local\_state
\]

严格禁止不同 sample 的细胞混合成同一 pool。跨样本比较发生在 sample-level aggregation 或 mixed model 中。

### Pool size

| 模式 | pool size | 用途 | 风险 |
|---|---:|---|---|
| Fast | 100–300 | 快速筛选 | 稀释稀有状态 |
| Balanced | 30–100 | 默认开发/发表 | 需要足够细胞数 |
| High-res | 10–30 | 稀有状态探索 | dropout 和 QP 不稳定 |

### Pool 表达

\[
X^{RNA}_{g,p}=\frac{1}{|C_p|}\sum_{c\in C_p}Z_{RNA}(g,c)
\]

\[
X^{ATAC}_{e,p}=\frac{1}{|C_p|}\sum_{c\in C_p}A(e,c)
\]

\[
M_{f,p}=\frac{1}{|C_p|}\sum_{c\in C_p}M_{f,c}
\]

### Pool 质量阈值

建议报告而非一刀切过滤：

```text
n_cells_per_pool
RNA_depth_per_pool
ATAC_depth_per_pool
metabolic_gene_detection_rate
GPR_gene_detection_rate
low_power_pool_flag
```

建议参考阈值：

```text
metabolic_gene_detection_rate ≥ 0.30
pathway-specific GPR gene detection_rate ≥ 0.50
```

不足时标记 low_power，而不是跨样本合并。

---

## Module 2：Gene-level multiome support

### 2.1 Enhancer support

\[
R_{ATAC}(g,p)=\sum_{e\in Enh(g)}w_{e,g}X^{ATAC}_{e,p}
\]

其中 \(w_{e,g}\) 来自 peak-gene link。必须保留：

```text
link_source
link_weight
link_confidence
distance_to_TSS
```

### 2.2 TF/motif support

\[
TFAct(f,p)=rank\_avg(RNA(f,p),MotifDev(f,p))
\]

\[
R_{TF}(g,p)=aggregate_{f\in TF(g)} TFAct(f,p)
\]

### 2.3 Gene-wise percentile concordance

不能在同一 pool 内比较不同 gene rank，因为会混入基因间 baseline 差异。应对每个 gene 跨 pool 计算 percentile：

\[
P_{RNA}(g,p)=percentile_p\left(X^{RNA}_{g,p}\right)
\]

\[
P_{ATAC}(g,p)=percentile_p\left(R_{ATAC}(g,p)\right)
\]

\[
Concord_{RA}(g,p)=1-|P_{RNA}(g,p)-P_{ATAC}(g,p)|
\]

同理：

\[
Concord_{RT}(g,p)=1-|P_{RNA}(g,p)-P_{TF}(g,p)|
\]

### 2.4 Fisher z shrinkage reliability

\[
\rho_{RA}(g)=Spearman(X^{RNA}_{g,\cdot},R_{ATAC}(g,\cdot))
\]

先 clip：

\[
\rho^{clip}=\min(0.999,\max(-0.999,\rho))
\]

Fisher z：

\[
z_g=arctanh(\rho^{clip}_g)
\]

收缩系数：

\[
\lambda_g=\frac{n_{pool}-3}{n_{pool}-3+n_0}
\]

默认：

\[
n_0=30
\]

\[
\rho^{shrink}_g=tanh(\lambda_gz_g)
\]

正相关作为支持：

\[
Rel_{RA}(g)=\max(0,\rho^{shrink}_{RA}(g))
\]

负相关保留为 discordance：

\[
Discord_{RA}(g)=|\min(0,\rho^{shrink}_{RA}(g))|
\]

### 2.5 Gene confidence

\[
Conf(g,p)=0.30\cdot Concord_{RA}(g,p)Rel_{RA}(g)
+0.20\cdot Concord_{RT}(g,p)Rel_{RT}(g)
+0.20\cdot Det_{RNA}(g,p)
+0.15\cdot LinkConf(g,p)
+0.15\cdot QC(p)
\]

约束：

\[
Conf(g,p)\in[0,1]
\]

**解释限制**：该 confidence 是 heuristic support score，不是概率。

---

## Module 3：GPR-aware reaction capacity potential

### 3.1 Gene score

对每个 gene 在 pool 间做 robust z-score：

\[
z_{g,p}=\frac{X^{RNA}_{g,p}-median_p(X^{RNA}_{g,p})}{MAD_p(X^{RNA}_{g,p})+\epsilon}
\]

转为 0–1 score：

\[
s_{g,p}=\sigma(z_{g,p})=\frac{1}{1+e^{-z_{g,p}}}
\]

### 3.2 Promiscuity correction

默认：

\[
s'_{g,p}=\frac{s_{g,p}}{\sqrt{N_{rxn}(g)}}
\]

必须做 sensitivity：

```text
none:       s' = s
sqrt mode:  s' = s / sqrt(N_rxn)
linear:     s' = s / N_rxn   # METAFlux-like
```

选择 \(1/\sqrt{N}\) 的理由：比 METAFlux 的 \(1/N\) 更温和，减少 hub gene 被过度惩罚；但必须报告 sensitivity。

### 3.3 Dropout-aware correction

使用 cell-type/substate-specific baseline，而不是 global baseline：

\[
m_{g,c}=median_{p\in celltype/substate\ c}(s'_{g,p})
\]

\[
\tilde{s}_{g,p}=q_{g,p}s'_{g,p}+(1-q_{g,p})m_{g,c(p)}
\]

其中 \(q_{g,p}\) 是 pool 内 gene detection rate。

若使用 global baseline，真实低表达细胞类型会被高表达细胞类型拉高，导致 capacity 膨胀。因此 global baseline 只能作为 sensitivity。

### 3.4 GPR parsing

所有 GPR 解析为 OR-of-AND groups：

```text
(g1 AND g2) OR g3 OR (g4 AND g5)
```

转为：

```text
AND_group_1 = {g1, g2}
AND_group_2 = {g3}
AND_group_3 = {g4, g5}
```

### 3.5 AND aggregation：Boltzmann-weighted average

不要称为标准 LogSumExp softmin。推荐名称：

```text
Boltzmann-weighted average
soft-minimum via weighted average
```

公式：

\[
w_{g,p}=\frac{\exp(-\tilde{s}_{g,p}/\tau)}{\sum_{h\in G_{r,k}}\exp(-\tilde{s}_{h,p}/\tau)}
\]

\[
C_k(r,p)=\sum_{g\in G_{r,k}}w_{g,p}\tilde{s}_{g,p}
\]

性质：

\[
\min_g \tilde{s}_{g,p}\le C_k(r,p)\le mean_g\tilde{s}_{g,p}
\]

因此它比 hard min 更抗 dropout，但会高于最小亚基，不是严格 bottleneck。

### 3.6 τ 参数验证

若希望两个亚基 score 差值为 \(\Delta=0.2\) 时，低 score 亚基的权重为高 score 亚基的 \(R=12\) 倍：

\[
\frac{w_{low}}{w_{high}}=\exp(\Delta/\tau)=R
\]

\[
\tau=\frac{\Delta}{\ln R}=\frac{0.2}{\ln 12}=0.0805
\]

因此默认：

\[
\tau=0.08
\]

必须报告 sensitivity grid：

\[
\tau\in\{0.03,0.05,0.08,0.12,0.16,0.20\}
\]

但 sensitivity 只在 Layer 1 capacity 上运行，不重复全量 QP。

### 3.7 OR aggregation

\[
C_{raw}(r,p)=\sum_{k=1}^{K_r}C_k(r,p)
\]

OR additive raw 保留 isoenzyme 累加效应，避免 Noisy-OR 或 saturation sum 对单 isoenzyme 产生不自然压缩。

### 3.8 Q95 relative calibration

\[
C_{rel}(r,p)=\min\left(1,\frac{C_{raw}(r,p)}{Q_r+\epsilon}\right)
\]

Q 值选择：

\[
Q_r=
\begin{cases}
Q_{0.95,r}^{direct}, & n_{pool}\ge100 \\
\rho Q_{0.95,r}^{stratum}+(1-\rho)Q_{0.95}^{global}, & 60\le n_{pool}<100 \\
Q_{0.90,r}^{stratum}, & 20\le n_{pool}<60 \\
Q_{0.95}^{global}, & n_{pool}<20
\end{cases}
\]

\[
\rho=\frac{n_{pool}}{n_{pool}+n_0}, \quad n_0=80
\]

即使 \(n_{pool}\ge100\)，Q95 仍由最高 5% observations 决定，因此必须输出 bootstrap CI：

```text
q95_bootstrap_CI_low
q95_bootstrap_CI_high
q95_CI_width
q95_unstable_flag = CI_width > 0.20
```

### 3.9 Absolute low flag

\[
low_r=I[median_p(C_{raw}(r,p))<\theta_{abs}]
\]

默认：

\[
\theta_{abs}=0.05
\]

该 flag 只能作为 diagnostic，不作为硬过滤。原因是 \(C_{raw}\) 受 isoenzyme 数、GPR group 数、promiscuity correction 和 cell type baseline 影响，不能跨 reaction 解释为绝对酶容量。

### 3.10 Reaction confidence

\[
RConf(r,p)=weighted\_median\{Conf(g,p):g\in GPR(r)\}
\]

无 GPR reaction：

```text
C_rel = 1.0 for network permissiveness
RConf = 0
flag = non_GPR
```

Transport reaction：

```text
flag = transport
RConf *= 0.5 if GPR incomplete
prioritize for Layer 3 exact scoring
```

---

## Module 4：GEM preprocessing 与 network integrity checks

### 4.1 GEM 选择

默认使用 Human-GEM / Human1，并记录：

```text
GEM_name
GEM_version
number_of_reactions
number_of_metabolites
number_of_genes
number_of_compartments
GPR_coverage
transport_reaction_fraction
exchange_reaction_fraction
```

Human-GEM 优点：GPR 体系较完整，Human1/METAFlux 使用过，stoichiometric consistency 较好。缺点：transport reaction 多、GPR 不完整、区室化复杂。

### 4.2 Flux consistency check

在 QP 前必须检查：

```text
blocked_reactions
dead_end_metabolites
mass_charge_imbalance_if_available
flux_consistent_fraction
```

Blocked reactions 不应进入 Layer 3 demand scoring，除非用户强制指定并明确标记。

### 4.3 Reversibility 与 v+/v- 分解

对可逆反应：

\[
v_r=v_r^+-v_r^-,\quad v_r^+\ge0,\quad v_r^-\ge0
\]

对不可逆反应：

\[
v_r^- = 0
\]

必须遵循 GEM 中的 reversibility / lower-bound / upper-bound 标注。

### 4.4 Futile cycle 诊断

即使有 \(L_2\) 正则化，也应记录：

```text
futile_cycle_flag[r,p] = (v_plus > eps and v_minus > eps)
futile_cycle_fraction_per_pool
```

### 4.5 区室化

必须保留 compartment 信息：

```text
reaction_compartment
metabolite_compartment
transport_flag
mitochondrial_flag
peroxisomal_flag
boundary_flag
```

Pathway aggregation 时，建议保留 compartment-specific reaction，而不是简单合并 cytosol/mitochondria 同名反应。

### 4.6 热力学与 loop diagnostics

默认不强制完整 thermodynamic FBA，因为缺少 metabolite concentration 与 ΔG 数据。但必须提供：

```text
loop_reaction_flag
loopless_FVA_option
thermodynamic_warning
```

若用户提供 ΔG 或 concentration range，可启用 thermodynamic constraint mode。

---

## Module 5：Network-constrained QP feasibility

### 5.1 推荐默认：Mode A soft-expression penalty

为了避免 double penalty，默认不把 RNA-derived capacity 直接作为 hard upper bound。GEM/nutrient 决定基本 bounds，expression-derived capacity 进入 objective penalty。

变量：

\[
v_r=v_r^+-v_r^-
\]

质量平衡：

\[
S(v^+-v^-)=0
\]

bounds：

\[
0\le v_r^+\le U^{model,+}_{r,p}
\]

\[
0\le v_r^-\le U^{model,-}_{r,p}
\]

penalty：

\[
penalty_{r,p}=-\log(C_{rel}(r,p)+10^{-6})
\]

baseline QP：

\[
\min_{v^+,v^-}\sum_r penalty_{r,p}(v_r^++v_r^-)+\lambda\left(\|v^+\|_2^2+\|v^-\|_2^2\right)
\]

subject to：

\[
S(v^+-v^-)=0
\]

\[
v_{ATPM}\ge \theta_{ATPM}
\]

\[
l^{model}\le v\le u^{model}
\]

推荐：

\[
\lambda\in[10^{-4},10^{-3}]
\]

### 5.2 Sensitivity：Mode B hard capacity bound

仅作为 sensitivity：

\[
U_{r,p}=U_{max}\cdot C_{rel}(r,p)
\]

Mode B 不应作为主分析，因为低 \(C_{rel}\) reaction 会同时被 hard bound 限制和 penalty 惩罚。

### 5.3 Exchange modes

至少运行：

| 模式 | 设置 | 解释 |
|---|---|---|
| conservative | 只开放核心营养、氧气、无机盐等 | 无 nutrient profile 默认 |
| medium-informed | 使用实验培养基或组织营养 profile | 推荐主分析 |
| liberal-penalized | 开放更多 exchange，但 uptake 加 penalty | sensitivity |

无空间组学、培养基、代谢物交换实测时，不建模 TME community exchange。

### 5.4 Biomass mode

仅用于明确增殖细胞或恶性细胞 sensitivity：

\[
\min \sum_r penalty_{r,p}|v_r|+\lambda\|v\|_2^2-\alpha v_{biomass}
\]

\[
\alpha=10^3\sim10^4
\]

biomass mode 与 ATP maintenance mode 必须分开输出，不能混合解释。

---

## Module 6：Selected reaction-demand exact feasibility

### 6.1 Reaction selection

Layer 3 仅对 selected reactions 做 exact demand QP：

\[
Selected=R_{exchange}\cup R_{transport}\cup R_{topLayer1Diff}\cup R_{bottleneck}\cup R_{highConf}\cup R_{user}\cup R_{validation}
\]

默认：

```text
N_selected = 300–800 for Balanced mode
N_selected ≤ 1500 for Accurate mode
```

必须报告：

\[
fraction_{L3}=\frac{N_{selected}}{N_{total}}
\]

### 6.2 Demand constraint

\[
\delta_{r,p}=\max(\delta_{min},0.01\cdot U^{model}_{r,p})
\]

不可逆 reaction：

\[
v_r^+\ge \delta_{r,p}
\]

可逆 reaction：分别测试 forward 和 reverse：

\[
v_r^+\ge \delta_{r,p}
\]

\[
v_r^-\ge \delta_{r,p}
\]

取较低 normalized cost，并记录：

```text
preferred_direction
forward_cost_norm
reverse_cost_norm
direction_ambiguous_flag
```

### 6.3 Normalized feasibility

\[
\Delta Obj_{r,p}=Obj_{demand}(r,p)-Obj_{baseline}(p)
\]

先截断数值误差：

\[
\Delta Obj^+_{r,p}=\max(0,\Delta Obj_{r,p})
\]

归一化：

\[
CostNorm_{r,p}=\frac{\Delta Obj^+_{r,p}}{\delta_{r,p}+\epsilon}
\]

\[
Feasibility_{norm}(r,p)=\frac{1}{1+CostNorm_{r,p}}
\]

同时保留 raw：

```text
DeltaObj_raw
Feasibility_raw = 1/(1+DeltaObj_positive)
Feasibility_norm
```

跨 reaction 排名使用 `Feasibility_norm`，同一 reaction 跨 pool 比较可同时参考 raw 与 norm。

---

## Module 7：Limited FVA

### 7.1 ATP maintenance mode

对 selected / uncertain reactions 运行：

\[
\min / \max \quad v_r
\]

subject to：

\[
S v=0
\]

\[
l\le v\le u
\]

\[
Obj(v)\le Obj^*+\eta_{abs}
\]

或当 \(Obj^*>0\) 且不接近 0：

\[
Obj(v)\le (1+\eta)Obj^*
\]

推荐：

```text
eta_abs = 1e-4 to 1e-3 after objective scaling
eta = 0.05 to 0.10
```

### 7.2 Biomass mode

不能使用：

\[
Obj(v)\le(1+\eta)Obj^*
\]

因为 \(Obj^*\) 可能为负。

应先求最大 biomass：

\[
v^*_{biomass}=\max v_{biomass}
\]

FVA 时约束：

\[
v_{biomass}\ge(1-\eta)v^*_{biomass}
\]

再对 selected reaction 做 min/max。

输出：

```text
FVA_min
FVA_max
FVA_width
direction_uncertain_flag
loop_flag
```

---

## Module 8：Pathway 与 meta-reaction aggregation

### 8.1 近常数过滤

\[
range_r=\max_p Score(r,p)-\min_p Score(r,p)
\]

若：

\[
range_r<10^{-3}
\]

标记为：

```text
non_informative
```

### 8.2 Pathway score

\[
PathwayScore(k,p)=weighted\_median\{Feasibility_{norm}(r,p):r\in Pathway(k)\}
\]

权重：

\[
w_{r,p}=RConf(r,p)\cdot I(not\ low\_r)
\]

低功效标记：

\[
lowPower(k,p)=I[n_{valid}(k,p)<3]
\]

### 8.3 Meta-reaction module

对 informative selected reactions 构建相关矩阵：

\[
D(r_i,r_j)=1-cor(Feasibility_{norm}(r_i,\cdot),Feasibility_{norm}(r_j,\cdot))
\]

使用 hierarchical clustering / Leiden on reaction graph 形成 meta-reaction modules。

---

## Module 9：Sample-aware statistics

### 9.1 默认：sample-level aggregation

\[
Y_{s,t,c,r}=median_{p\in(s,t,c)} Score(p,r)
\]

模型：

\[
Y_{s,t,c,r}=\beta_0+\beta_1Condition_s+\beta_2Batch_s+\beta_3CNV_s+\beta_4CellCycle_s+\epsilon_s
\]

### 9.2 Mixed model

当样本数足够且每 sample 有多个 pool：

\[
Score_{p,r}=\beta_0+\beta_1Condition_p+\beta_2PoolSize_p+\beta_3RNAdepth_p+\beta_4ATACdepth_p+u_{sample[p]}+\epsilon_p
\]

\[
u_{sample}\sim N(0,\sigma_s^2)
\]

最低要求：

```text
n_samples_per_group ≥ 3 for exploratory
n_samples_per_group ≥ 5–6 for more stable mixed model
```

### 9.3 多重检验

在每个 cell type 内分别进行 BH-FDR：

\[
q_r=BH(p_r)
\]

输出字段：

```text
reaction_id
cell_type
contrast
effect_size
p_value
FDR
n_samples_per_group
n_pools_per_group
scoring_layer
q95_mode
globally_low_flag
tau_sensitive_flag
exchange_sensitive_flag
FVA_width
pseudo_replication_warning
```

---

## Module 10：Signaling-anchored metabolic rewiring 可选模块

用于 EGFR、ERBB2、MYC、HIF1A、mTOR、TGFβ 等 signaling-metabolism 问题。

### 10.1 Signaling activity

不要只用单基因表达。以 EGFR 为例：

\[
EGFRActivity(p)=w_1 EGFR_{RNA}(p)+w_2 ERBB\ pathway(p)+w_3 MAPK/PI3K\ downstream(p)+w_4 AP1/ETS/MYC\ motif(p)
\]

默认建议使用 pathway/regulon score，而非 EGFR mRNA 单基因。

### 10.2 Association model

\[
Score(r,p)\sim EGFRActivity(p)+epithelial\_substate(p)+CNV(p)+CellCycle(p)+RNAdepth(p)+ATACdepth(p)+(1|sample)
\]

或 within-sample correlation + cross-sample meta-analysis：

\[
\rho_s(r)=cor(EGFRActivity_{p\in s},Score(r,p))
\]

\[
z_s=arctanh(\rho_s)
\]

跨样本检验 \(z_s\) 是否稳定偏离 0。

### 10.3 输出解释

可以说：

```text
EGFR-associated metabolic capacity / feasibility program
EGFR-associated candidate regulatory mediator
sample-stable association
```

不能说：

```text
EGFR causally drives this flux
TF is proven driver
reaction flux is activated by EGFR
```

---

## Module 11：Regulator candidate prioritization

### 11.1 输入

```text
TF RNA
motif deviation
enhancer accessibility
peak-gene links
target metabolic genes
reaction capacity / feasibility
sample metadata
```

### 11.2 Direct association

\[
Assoc_{direct}(f,r)=cor(TFAct(f,\cdot),Score(r,\cdot))
\]

### 11.3 Expression-adjusted association

先残差化 reaction score：

\[
Score(r,p)\sim TargetGeneExpr(r,p)+Batch_p+PoolSize_p+Sample_p
\]

取残差：

\[
ResidScore(r,p)
\]

再计算：

\[
Assoc_{adj}(f,r)=cor(TFAct(f,\cdot),ResidScore(r,\cdot))
\]

### 11.4 Collinearity control

对证据源做相关过滤：

\[
|cor(E_i,E_j)|>0.9
\]

保留更稳定或更上游的证据源。

### 11.5 Ranking

使用 robust rank aggregation：

```text
rank_direct_association
rank_adjusted_association
rank_motif_support
rank_enhancer_support
rank_sample_bootstrap_stability
rank_leave_one_sample_out_stability
```

等级：

| 等级 | 含义 |
|---|---|
| correlation-only | 只有相关 |
| motif-supported | 有 motif 支持 |
| enhancer-supported | 有 enhancer-gene link |
| multiome-supported | RNA+ATAC+motif 一致 |
| sample-stable | leave-one-sample-out 稳定 |
| perturbation-supported | 外部扰动验证支持 |

---

## Module 12：Cell-wise projection

默认不做逐细胞 QP。

### Pool assignment

\[
CellScore(c,r)=Score(pool(c),r)
\]

### WNN barycentric interpolation

\[
CellScore(c,r)=\sum_{k\in KNN(c)}w_{c,k}Score(pool_k,r)
\]

\[
\sum_k w_{c,k}=1
\]

该分数仅用于 UMAP 可视化，不能用于 p-value。

---

## 6. 计算验证细节与单元测试

### 6.1 τ 验证测试

目标：验证默认 \(\tau=0.08\) 的权重比逻辑。

输入：

```text
s_low = 0.4
s_high = 0.6
Delta = 0.2
R_target = 12
```

计算：

\[
R=\exp((s_{high}-s_{low})/\tau)
\]

当 \(\tau=0.08\)：

\[
R=\exp(0.2/0.08)=12.18
\]

通过条件：

```text
abs(R - 12) / 12 < 0.05
```

### 6.2 Boltzmann-weighted average 范围测试

对随机 scores：

\[
C_k=\sum_gw_gs_g
\]

通过条件：

```text
min(scores) <= C_k <= mean(scores)
```

如果使用 LogSumExp softmin：

\[
-\tau\log\sum_g e^{-s_g/\tau}
\]

该值可能小于 min，甚至为负，不用于主公式。

### 6.3 Fisher z shrinkage 测试

输入：

```text
rho = 1.0, -1.0, 0.8, -0.4
```

必须先 clip 到：

```text
[-0.999, 0.999]
```

通过条件：

```text
all finite(arctanh(rho_clip))
abs(rho_shrink) <= abs(rho_clip)
```

### 6.4 Dropout baseline 膨胀测试

模拟：

```text
celltype_A true s = 0.1
celltype_B true s = 0.8
global_median = 0.45
q = 0.3
```

global baseline：

\[
\tilde{s}=0.3\cdot0.1+0.7\cdot0.45=0.345
\]

cell-type baseline：

\[
\tilde{s}=0.3\cdot0.1+0.7\cdot0.1=0.1
\]

通过结论：global baseline 会把低表达状态膨胀 3.45 倍，因此主实现使用 cell-type/substate baseline。

### 6.5 Q95 稳定性测试

对每个 reaction bootstrap pools B=200：

```text
Q95_boot[b] = quantile(C_raw_boot, 0.95)
CI_width = Q97.5 - Q2.5
q95_unstable_flag = CI_width > 0.20
```

若 unstable，降级为 shrinkage/global Q。

### 6.6 双重惩罚测试

同一 toy network 比较：

```text
Mode A: soft penalty only
Mode B: hard bound + penalty
Mode C: hard bound only
```

预期：

```text
Mode B <= Mode A in feasible flux space
Mode B and C may be similar when low-capacity reaction is hard bottleneck
```

主分析使用 Mode A，Mode B 作 sensitivity。

### 6.7 Demand normalization 测试

比较：

```text
r1: delta=1,  DeltaObj=5
r2: delta=10, DeltaObj=5
```

raw feasibility 相同：

\[
1/(1+5)=0.167
\]

normalized feasibility：

\[
r1=1/(1+5/1)=0.167
\]

\[
r2=1/(1+5/10)=0.667
\]

结论：normalized feasibility 更符合单位 demand 成本逻辑。

### 6.8 FVA biomass constraint 测试

若：

\[
Obj^* = penalty - \alpha biomass = -498
\]

则：

\[
(1+0.05)Obj^*=-522.9
\]

该约束要求比最优解更优，数学错误。biomass mode 改用：

\[
v_{biomass}\ge(1-\eta)v^*_{biomass}
\]

---

## 7. 实现架构

建议主实现 Python，兼容 AnnData/MuData；R 作为接口层。

```text
regcompass/
  io/
    read_multiome.py
    read_gem.py
    write_outputs.py
  preprocessing/
    qc.py
    normalization.py
    pooling.py
    pseudobulk.py
  regulatory/
    peak_gene.py
    motif.py
    confidence.py
  gpr/
    parser.py
    aggregate.py
    calibration.py
    diagnostics.py
  gem/
    model.py
    consistency.py
    compartments.py
    reversibility.py
  qp/
    build.py
    baseline.py
    demand.py
    fva.py
    sensitivity.py
  stats/
    sample_level.py
    mixed_model.py
    permutation.py
    fdr.py
  regulators/
    association.py
    residualization.py
    rank_aggregation.py
  signaling/
    activity.py
    association.py
  visualization/
    project_cells.py
  report/
    diagnostics.py
    markdown.py
  tests/
    test_softmin.py
    test_q95.py
    test_fisher.py
    test_qp.py
    test_fva.py
```

### 7.1 核心矩阵

```text
gene_score: genes × pools
gene_confidence: genes × pools
reaction_capacity_L1: reactions × pools
reaction_confidence: reactions × pools
baseline_qp_diagnostics: pools × diagnostics
reaction_feasibility_L3: selected_reactions × pools
pathway_score: pathways × pools
regulator_ranking: TFs × reactions/pathways
cell_projection: cells × selected_reactions/pathways
```

### 7.2 Solver

推荐：

```text
OSQP for convex QP
HiGHS/scipy.optimize for LP/FVA
Gurobi optional for large-scale exact solve
```

OSQP 适合原因：

```text
S matrix sparse and stable across pools
only q/l/u bounds change
warm start possible
factorization/matrix structure reusable
infeasibility status available
```

### 7.3 加速策略

必须实现：

```text
sparse S matrix
reaction split cached once
QP matrix cache
warm start from baseline
parallel by pool or reaction
selected reactions only
skip blocked reactions
checkpoint every N solves
resume from checkpoint
solver status logging
```

### 7.4 计算成本估算

运行前必须输出：

\[
N_{QP}=N_{pool}\times(1+N_{selected})
\]

例如：

```text
1000 pools × (1 baseline + 500 demand) = 501000 QP solves
```

工具必须在 dry-run 中报告：

```text
estimated_QP_count
estimated_memory
estimated_runtime_range
recommended_parallel_jobs
```

---

## 8. 运行模式

| 参数 | Fast | Balanced | Accurate |
|---|---:|---:|---:|
| pool size | 100–300 | 30–100 | 10–50 |
| Layer 0 | 是 | 是 | 是 |
| Layer 1 | 是 | 是 | 是 |
| baseline QP | 否 | 是 | 是 |
| selected demand QP | 否 | 300–800 | ≤1500 |
| FVA | 否 | limited | expanded limited |
| τ sensitivity | Layer 1 only | Layer 1 only | Layer 1 + selected |
| exchange sensitivity | optional | required | required |
| statistics | descriptive | sample-level | sample-level + mixed/permutation |
| regulator ranking | rough | full | full + bootstrap |
| target | exploration | default publication | validation-focused |

---

## 9. 验证体系

### 9.1 内部负对照

```text
random GPR mapping
random peak-gene links matched by distance/accessibility
motif label permutation
condition label permutation within valid block
sample label permutation
expression-matched random gene sets
```

### 9.2 稳定性分析

```text
pool bootstrap
leave-one-sample-out
softmin τ grid
Q95 bootstrap CI
promiscuity mode sensitivity
exchange mode sensitivity
N_selected sweep
GPR perturbation 5% / 10% / 20%
Mode A vs Mode B QP comparison
```

### 9.3 外部验证优先级

1. matched metabolomics
2. 13C isotope tracing
3. Seahorse OCR/ECAR
4. targeted metabolic enzyme perturbation
5. TF perturbation with RNA+ATAC readout
6. proteomics / phosphoproteomics
7. COMPASS / METAFlux concordance

### 9.4 与参考工具的验证

- 与 scMetabolism：pathway rank correlation，仅作为 sanity check。
- 与 scFEA：module-level concordance，不要求 reaction-level 一致。
- 与 COMPASS：Layer 1 penalty / selected reaction direction consistency。
- 与 METAFlux：在相同 Human-GEM、nutrient profile 与 cell-type pseudobulk 上比较 selected reaction rank concordance。

---

## 10. Diagnostics report 必须包含

### Pool diagnostics

```text
n_cells_per_pool
n_pools_per_sample
n_pools_per_cell_type
low_power_pool_fraction
metabolic_gene_detection_rate
GPR_gene_detection_rate
```

### Multiome diagnostics

```text
RNA depth distribution
ATAC depth distribution
TSS enrichment
FRiP
RNA-ATAC concordance distribution
motif support distribution
negative_concordance_fraction
```

### GPR diagnostics

```text
GPR coverage
non_GPR_fraction
transport_reaction_fraction
multi_subunit_fraction
isoenzyme_fraction
promiscuous_gene_fraction
q95_mode_distribution
q95_unstable_fraction
globally_low_capacity_fraction
```

### GEM/QP diagnostics

```text
GEM version
blocked reaction fraction
flux consistent reaction fraction
baseline_infeasible_rate
demand_infeasible_rate
solver_status_distribution
median_QP_runtime
futile_cycle_fraction
exchange_sensitivity_score
fraction_L3_exact
```

### Statistics diagnostics

```text
n_samples_per_group
n_pools_per_group
sample_imbalance_warning
batch_condition_confounding_warning
pseudo_replication_warning
```

---

## 11. 默认参数表

| 参数 | 默认值 | 必须报告 | Sensitivity |
|---|---:|---:|---:|
| pool size | 30–100 | 是 | 10–30, 100–300 |
| AND aggregation | Boltzmann-weighted average | 是 | hard min, mean |
| τ | 0.08 | 是 | 0.03–0.20 |
| OR aggregation | additive raw | 是 | max, saturation sum |
| promiscuity | \(1/\sqrt{N_{rxn}}\) | 是 | none, \(1/N\) |
| dropout baseline | cell-type/substate median | 是 | global median |
| Q95 direct | n≥100 | 是 | bootstrap CI |
| Q95 shrinkage n0 | 80 | 是 | 40, 120 |
| Fisher shrinkage n0 | 30 | 是 | 10, 50 |
| absolute low θ | 0.05 flag only | 是 | 0.02, 0.10 |
| QP mode | soft penalty | 是 | hard bound sensitivity |
| λ L2 | 1e-4–1e-3 | 是 | 1e-5–1e-2 |
| baseline task | ATP maintenance | 是 | biomass mode |
| selected reactions | 300–800 | 是 | ≤1500 |
| demand δ | max(δmin, 0.01Umodel) | 是 | 0.005, 0.02 |
| FVA η | 0.05–0.10 | 是 | 0.01, 0.20 |
| KNN diffusion | off | 是 | trajectory-only sensitivity |
| statistics | sample-level | 是 | mixed model, blocked permutation |

---

## 12. 最终伪代码

```text
Input:
  RNA_counts, ATAC_counts, metadata,
  WNN_embedding, peak_gene_links, motif_deviation,
  Human-GEM, nutrient_profile(optional)

Step 0: QC and normalization
  QC RNA and ATAC
  normalize RNA with sample-aware residuals
  normalize ATAC with TF-IDF/LSI
  build or import WNN embedding

Step 1: Sample-aware micropooling
  for each sample × condition × cell_type:
    define local states in WNN space
    construct micropools without crossing sample
    compute RNA/ATAC/motif pseudobulk
    record pool diagnostics

Step 2: Gene-level support
  for each gene and pool:
    compute RNA score
    compute enhancer ATAC support
    compute TF/motif support
    compute gene-wise percentile concordance
    compute Fisher-z shrinkage reliability
    output gene confidence and discordance

Step 3: GPR capacity
  parse GPR into OR-of-AND groups
  compute robust z-score and sigmoid gene score
  apply promiscuity correction
  apply cell-type-specific dropout correction
  compute AND by Boltzmann-weighted average
  compute OR by additive raw sum
  apply Q95 calibration with bootstrap CI
  output C_raw, C_rel, RConf, q95 flags, low flags

Step 4: GEM preprocessing
  load Human-GEM
  check blocked reactions and flux consistency
  split reversible reactions respecting model bounds
  annotate compartments and transport reactions

Step 5: Baseline QP
  for each pool:
    construct model/nutrient bounds
    use C_rel as soft penalty
    solve ATP-maintenance baseline QP
    record objective, infeasibility, runtime, futile cycles

Step 6: Selected demand QP
  select exchange, transport, top Layer1, bottleneck, high-conf, user reactions
  for each selected reaction and pool:
    impose demand constraint
    solve demand QP with warm start
    compute DeltaObj, CostNorm, Feasibility_norm
    record preferred direction and infeasible flag

Step 7: Limited FVA
  run FVA for selected or uncertain reactions
  use correct objective/biomass constraints
  output width and direction uncertainty

Step 8: Aggregation
  filter near-constant reactions
  compute pathway weighted median score
  cluster meta-reaction modules

Step 9: Statistics
  aggregate pools to sample level
  run sample-level model or mixed model
  adjust BH-FDR
  output pseudo-replication diagnostics

Step 10: Regulator ranking
  compute direct and expression-adjusted TF-reaction associations
  control collinearity
  rank by RRA and sample stability

Step 11: Optional signaling module
  compute EGFR/ERBB/MYC/HIF/mTOR activity
  test signaling-associated metabolic rewiring
  report candidate mediators only

Step 12: Cell projection
  assign pool scores to cells or WNN barycentric interpolation
  mark as visualization-only

Output:
  all matrices, diagnostics report, sensitivity summary, regulator ranking, projected scores
```

---

## 13. 最终结论

报告的主要判断是正确的，但需要将若干数学和工程细节进一步收紧。调整后的 RegCompass-Multiome 应被实现为：

```text
multiome-supported
GPR-aware
sample-aware
reaction-capacity-oriented
selected-QP-feasibility-oriented
diagnostic-rich
```

它的核心贡献不是“更准确预测真实 flux”，而是：

> 在 paired scRNA+scATAC multiome 中，将 RNA-derived GPR reaction capacity、ATAC/motif regulatory support、Human-GEM network constraints 和 sample-aware statistics 组合起来，识别具有多组学支持和网络可行性的代谢反应/通路/候选调控因子。

它必须持续声明以下边界：

```text
capacity potential ≠ enzyme activity
feasibility ≠ true flux
ATAC support ≠ reaction activation
regulator ranking ≠ causal driver
pool-level score ≠ independent single-cell flux
selected exact QP ≠ full reaction coverage
```

若按本文规格开发，该工具具备可实现性和发表潜力；若忽略上述边界，则会产生严重的生物学过度解释和数学伪精确。

---

## 14. 参考资料

1. COMPASS documentation: https://yoseflab.github.io/Compass/
2. COMPASS Settings: https://yoseflab.github.io/Compass/Compass-Settings.html
3. Wagner A, Wang C, Fessler J, et al. Metabolic modeling of single Th17 cells reveals regulators of autoimmunity. Cell. 2021.
4. Huang Y, Mohanty V, Dede M, et al. Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux. Nature Communications. 2023.
5. METAFlux GitHub: https://github.com/KChen-lab/METAFlux
6. Alghamdi N, Chang W, Dang P, et al. A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data. Genome Research. 2021.
7. scFEA GitHub: https://github.com/changwn/scFEA
8. scMetabolism GitHub: https://github.com/wu-yc/scMetabolism
9. Robinson JL, Kocabaş P, Wang H, et al. An atlas of human metabolism. Science Signaling. 2020.
10. Human-GEM repository: https://github.com/SysBioChalmers/Human-GEM
11. OSQP documentation: https://osqp.org/docs/
