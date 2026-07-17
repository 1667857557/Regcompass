# RegCompassR PR #110 各部分修改总结

## 1. 文档范围

本文总结 PR **#110 — Correct biological evidence, GPR logic, and inference semantics** 对 RegCompassR 的主要修改。

- 仓库：`1667857557/Regcompass`
- 分支：`agent/fix-math-biological-architecture`
- 目标分支：`main`
- 版本：`1.4.0 → 1.4.1`
- PR 状态：Draft、可合并
- CI 状态：源码解析、依赖安装、包安装和完整 `testthat` 测试均通过
- 公共 API：继续保持 4 个入口不变

```r
rc_prepare_human2_gem()
rc_make_medium_scenarios()
rc_run_regcompass()
rc_run_regcompass_one_shot()
```

本次修改的目标不是把 RNA+ATAC 推断升级为真实定量通量，而是修复原流程中可能导致错误代谢支持解释的数学退化、GPR 逻辑、调控证据、培养基语义、统计单位和最终评分语义。

---

# 2. 十二个问题的最终处理状态

| 编号 | 原问题 | 最终核查 | 本次处理 |
|---|---|---|---|
| 1 | 全零或恒定基因变成高容量 | 确定存在 | 改为零保留的绝对支持函数；Q95 不再作为 LP 容量 |
| 2 | 绝对表达和相对状态混为一体 | 确定存在 | 分离绝对支持与相对状态输出 |
| 3 | GPR 复杂逻辑被丢弃 | 导入路径确定存在 | 增加递归 Boolean parser；解析失败直接终止 |
| 4 | 多亚基酶 soft-min 被其他亚基补偿 | 数学上存在 | canonical 默认改为 `AND = min` |
| 5 | 同工酶求和造成 isoenzyme-count bias | 确定存在 | canonical 默认改为 `OR = max` |
| 6 | promiscuity 权重任意 | 默认实际为 \(1/\sqrt K\) | canonical 默认改为 `none` |
| 7 | Pando confidence 实际主要是 peak accessibility | 基本成立 | 改为有符号 TF × peak × gene 调控支持 |
| 8 | 浓度被等同于摄取速率 | 原表述不准确 | 将 generic medium 改名为技术性 permissive baseline |
| 9 | shared-TF 投影丢失方向和符号 | 确定存在 | 保留 regulator、target、sign；冲突边不参与组件合并 |
| 10 | FASTCORE support 被赋予低 penalty | canonical 主流程原本不成立 | 禁止旧 helper 默认进行 structural support 覆盖 |
| 11 | metacell 被当作生物重复 | 主流程未直接检验，但易误用 | 默认推断单位改为 sample × cell type |
| 12 | MAD-sigmoid 变成伪概率 | 转换确定存在 | raw penalty 设为主输出；相对值改为经验排名 |

---

# 3. 表达证据与反应容量修正

## 3.1 原有数学退化

原流程首先对每个基因进行跨 metacell 的稳健标准化：

\[
z_{gm}
=
\frac{x_{gm}-\operatorname{median}(x_g)}
{\max\{\operatorname{MAD}(x_g),\operatorname{IQR}(x_g)/1.349,0.05\}}
\]

再计算：

\[
G_{gm}=\sigma(z_{gm})
\]

如果一个基因在所有 metacell 中均为 0：

\[
x_{gm}=0
\Rightarrow
z_{gm}=0
\Rightarrow
G_{gm}=0.5
\]

随后按每条反应自己的 Q95 归一化：

\[
C^{relative}_{rm}
=
\min\left(
\frac{C^{raw}_{rm}}
{Q_{0.95}(C^{raw}_{r})+\epsilon},
1
\right)
\]

对于恒定的 0.5：

\[
C^{relative}_{rm}
\approx
\frac{0.5}{0.5+\epsilon}
\approx1
\]

最终：

\[
-\log C^{relative}_{rm}\approx0
\]

即全零基因可能产生近乎最低的 LP penalty。

## 3.2 新的零保留支持函数

canonical 流程现在使用：

\[
A=\frac{x}{x+\kappa}
\]

其中：

- \(x\) 是已经输入流程的非负标准化信号；
- RNA 中 \(x=\log(1+\mathrm{CPM})\)；
- ATAC 中 \(x\) 是对应输入尺度上的非负归一化可及性；
- \(\kappa>0\) 是半饱和参数。

因此：

\[
x=0\Rightarrow A=0
\]

且：

\[
x>0\Rightarrow0<A<1
\]

该值的正确名称是：

> bounded molecular support score

它不是概率、酶浓度、\(V_{\max}\)、物理反应容量或真实代谢通量。

## 3.3 为什么直接使用输入信号而不是反变换为 CPM

初步方案曾考虑：

\[
A=
\frac{\exp(x)-1}
{\exp(x)-1+\kappa}
\]

但 metacell CPM 往往较高，该变换容易快速饱和；同时它不适用于 TF-IDF 或其他 ATAC 信号。

最终统一为：

\[
A=\frac{x}{x+\kappa}
\]

优点：

1. 零点严格保留；
2. RNA 与 ATAC 可分别在自身已归一化尺度上使用；
3. 避免对高 CPM 过早饱和；
4. 保持单调性；
5. 不虚构 RNA 与 ATAC 的共同物理单位。

## 3.4 绝对支持与相对状态分离

现在保留两套不同语义的结果。

### 绝对支持

用于 LP penalty：

\[
A^{abs}_{gm}
\]

\[
C^{abs}_{rm}
\]

### 相对状态

用于描述同一基因或同一反应在不同单位间的位置：

\[
R^{relative}_{gm}
\]

\[
C^{within-reaction}_{rm}
\]

Q95 归一化仍可输出为诊断项：

```text
C_within_reaction_relative
```

但不再用于 LP 主容量。

兼容字段 `C_rel` 目前保存的是 bounded absolute reaction support，而不是旧版 reaction-wise Q95 ratio。

---

# 4. GPR 逻辑修正

## 4.1 递归 Boolean parser

旧 parser 主要支持平坦规则，例如：

\[
(g_1\land g_2)\lor g_3
\]

但不能安全处理：

\[
g_1\land(g_2\lor g_3)
\]

新 parser 支持嵌套括号、`AND` 优先于 `OR`、任意合法组合，并转换为现有下游函数可使用的 DNF 表示。

例如：

\[
g_1\land(g_2\lor g_3)
\Rightarrow
(g_1\land g_2)\lor(g_1\land g_3)
\]

为避免组合爆炸，DNF 展开设有最大 term 数限制。超过限制时要求用户提供结构化 long-table GPR。

## 4.2 取消静默丢弃

旧 Human-GEM 导入可将 parser error 转换为空列表，从而静默丢弃反应的 GPR。

现在任何失败都会报告 reaction ID、原始 GPR 字符串和具体 parser 错误，并终止导入。

## 4.3 long-table GPR 的约束

对于 long-table 输入，必须显式提供：

```text
reaction_id
and_group_id
gene
```

缺少 `and_group_id` 时直接停止，避免把一个复合物中的多个必需亚基错误解释为多个独立同工酶。

---

# 5. 多亚基酶、同工酶和 promiscuity

## 5.1 多亚基复合物

canonical 默认从 Boltzmann soft-min 改为：

\[
C^{AND}_{rm}
=
\min_{g\in AND_r} A_{gm}
\]

生物学含义是完整酶复合物的支持由最弱的必需亚基限制。

## 5.2 同工酶

canonical 默认从：

\[
\frac{\sum_{k=1}^{K}C_k}{\sqrt K}
\]

改为：

\[
C^{OR}_{rm}
=
\max_k C_{km}
\]

原因是 RNA 不是 molar enzyme abundance，且缺少 \(k_{\mathrm{cat}}\)。求和会产生 isoenzyme-count bias。

如果未来获得蛋白丰度和动力学参数，物理容量应使用：

\[
V_{\max,r}
=
\sum_k k_{\mathrm{cat},k}E_k
\]

## 5.3 promiscuity

canonical 默认：

```r
promiscuity_mode = "none"
```

旧的 \(1/\sqrt K\) 和 \(1/K\) 仅保留为显式敏感性分析。

真正的酶资源竞争应表达为：

\[
\sum_{r\in R(g)}
\frac{|v_r|}
{k_{\mathrm{cat},gr}}
\le E_g
\]

---

# 6. Pando 调控证据修正

## 6.1 原有证据的局限

旧实现主要以：

\[
|\hat\beta_e|\sqrt{R_e^2}
\]

加权 peak accessibility，丢弃 coefficient 正负号，也未要求 TF 在当前单位中活跃。

## 6.2 新的有符号调控支持

对于 TF \(t\)、peak \(p\)、target gene \(g\) 的 edge \(e\)：

\[
u_{em}
=
\operatorname{sign}(\hat\beta_e)
T_{tm}
A_{pm}
\]

其中 \(T_{tm}\) 为 TF RNA 支持，\(A_{pm}\) 为 peak accessibility 支持。

基因层调控支持：

\[
F_{gm}
=
\operatorname{clip}_{[0,1]}
\left(
0.5+
0.5\sum_{e\rightarrow g}
w_eu_{em}
\right)
\]

解释：

- \(F=0.5\)：中性或缺少当前调控证据；
- \(F>0.5\)：正向调控支持；
- \(F<0.5\)：活跃抑制。

## 6.3 多组学 penalty 的新语义

RNA 仍是酶可用性的主要证据。缺失或中性调控不增加 penalty，活跃抑制增加非负 penalty，正向调控不会产生负 penalty。

---

# 7. shared-TF 投影修正

shared-TF 和 direct-regulatory 边现在保留 regulator、target、regulator set、coefficient sign、signed projection weight 和 concordant/discordant/mixed relation。

组件构建采用更保守规则：

- direct TF→target metabolic edge：可参与组件构建；
- concordant shared-TF edge：可参与组件构建；
- discordant 或 mixed shared-TF edge：仅保留诊断，不再合并 biological component。

完整机制解释仍应回到原始：

\[
TF\rightarrow peak\rightarrow gene\rightarrow reaction
\]

有向图。

---

# 8. Medium 语义修正

原 `blood_like` 和 `culture_like` 对全部 exchange reaction 开放摄取，使用统一 bound，但不读取真实浓度或 uptake/secretion rate。

canonical 默认改为：

```text
permissive_all_exchange
```

并标记：

```text
assumption_level = technical_upper_bound
```

旧名称仍可兼容，但会提示它们不是 curated biological medium。

新输出明确：

```text
concentration_used_for_rate_bound = FALSE
```

Control 与 Cre 原则上应使用相同 medium，除非存在独立环境测量支持 condition-specific exchange constraints。

---

# 9. FASTCORE structural support 修正

canonical 主流程原本没有因 `fastcore_support = TRUE` 就统一赋予低 penalty。

本次主要修复旧 helper 的误用风险：默认禁止 structural support 覆盖 biological evidence penalty。只有显式设置：

```r
allow_structural_support_override = TRUE
```

才允许旧行为。

原则为：

\[
\text{structural support}
\neq
\text{biological evidence}
\]

---

# 10. 推断单位修正

canonical 默认从：

```r
unit = "metacell"
```

改为：

```r
unit = "sample_celltype"
```

metacell 仍可用于样本内异质性和探索性可视化，但不能直接充当独立生物重复。

差异分析 helper 在 raw penalty 存在时优先使用 penalty，并在 sample × cell type 层聚合。

---

# 11. LP penalty 与最终 score 修正

两阶段方向性 LP 保持不变。

第一阶段：

\[
v_r^{max}
=
\max_v d_rv_r
\]

第二阶段：

\[
P^*_{rm}
=
\min_v
\sum_jp_{jm}|v_j|
\]

约束：

\[
Sv=0,\qquad
l\le v\le u,\qquad
d_rv_r\ge\omega v_r^{max}
\]

主输出解释为：

> 在指定 GEM、medium、方向和目标需求下，使目标反应可行所需的最小多组学证据不一致代价。

raw penalty 成为主推断输出。

旧 MAD-sigmoid 被替换为：

\[
R_{rm}
=
1-
\frac{\operatorname{rank}(P_{rm})-1}
{n_r-1}
\]

语义为：

```text
within_target_relative_penalty_rank_not_probability
```

恒定 target 返回 `NA` 并标记 `noninformative_target = TRUE`。

---

# 12. 导出、版本和 CI 工程修改

- `DESCRIPTION`：版本更新为 1.4.1；
- `NEWS.md`：新增修复说明；
- 新增架构说明文档；
- 保持公共 API 不变；
- gzip 输出连接改为显式关闭；
- CI 无论成功或失败都上传测试诊断 artifact。

最终 GitHub Actions run 159 中源码解析、依赖安装、包安装、完整 `testthat` 和 artifact 上传均通过。

两个依赖 Seurat 的既有 SuperCell membership 测试在最小 CI 环境中跳过；所有实际执行测试均通过。

---

# 13. 主要变更文件

## 核心实现

- `R/zzz_architecture_correctness.R`
- `R/zzzz_architecture_hotfixes.R`
- `R/zzzzz_signed_projection.R`
- `R/zzzzzz_signal_scale.R`
- `R/export.R`

## 文档和版本

- `DESCRIPTION`
- `NEWS.md`
- `docs/architecture-corrections.md`
- `.github/workflows/fastcore-checks.yaml`

## 测试

- `tests/testthat/test-architecture-correctness.R`
- `tests/testthat/test-signed-projection-components.R`
- `tests/testthat/test-medium-human-gem.R`
- `tests/testthat/test_gpr_parser.R`
- `tests/testthat/test_global_metacell_workflow.R`
- `tests/testthat/test_integration.R`
- `tests/testthat/test_microcompass_corrections.R`
- `tests/testthat/test_public_api.R`

---

# 14. 新增关键测试不变量

1. \(x=0\Rightarrow A=0\)；
2. 恒定正信号保留正支持；
3. nested GPR 正确展开；
4. `AND = min` 保留必需亚基瓶颈；
5. `OR = max` 避免 isoenzyme-count inflation；
6. missing regulation 中性，active repression 增加 penalty；
7. structural support 不得静默覆盖 biological penalty；
8. constant penalty 不生成伪相对 score；
9. permissive medium 明确标记为技术性假设；
10. shared-TF 保留 sign，并阻止冲突边合并组件；
11. sample × cell type 为默认推断单位。

---

# 15. 使用兼容性变化

## 15.1 新旧 score 不可直接混合

旧 score 基于 MAD-sigmoid；新版相对值基于 penalty rank。二者不能直接用于纵向合并比较。

## 15.2 `C_rel` 语义改变

为了减少调用链破坏，字段名保留，但其内容现在是 bounded absolute support。

旧的 reaction-wise ratio 位于：

```text
C_within_reaction_relative
```

## 15.3 Medium 名称变化

推荐：

```r
medium_scenario = "permissive_all_exchange"
```

旧 `blood_like` 仍可调用，但会提示其非 curated medium。

## 15.4 metacell 分析

需要显式：

```r
unit = "metacell"
```

并仅作为探索性结果解释。

---

# 16. 本次修改仍未解决的问题

1. RNA+ATAC 仍不能识别真实绝对通量；
2. 尚未加入 enzyme allocation constraints；
3. 尚未加入完整 thermodynamic/loopless constraints；
4. Pando 尚未完整加入 bootstrap edge frequency、cross-validation 和 sample random effects；
5. `permissive_all_exchange` 仍只是上界敏感性场景；
6. 小鼠数据最终仍应使用 Mouse-GEM 或经过严格验证的 orthology 转换模型。

真实定量 flux 仍需要 uptake/secretion measurements、extracellular metabolomics、蛋白或酶约束、\(k_{\mathrm{cat}}\) 和最好匹配的同位素示踪。

---

# 17. 修改后的正确方法学定位

修改后 RegCompassR 最准确的定位是：

> 基于 RNA、ATAC、Pando、GPR、GEM、medium 和稳态质量平衡约束，计算方向性反应或代谢任务的多组学兼容支持以及最小证据不一致代价。

核心关系：

\[
\text{RNA/ATAC support}
\rightarrow
\text{GPR reaction evidence}
\rightarrow
\text{evidence penalty}
\rightarrow
\text{directional constrained LP}
\rightarrow
\text{minimum evidence-discordance penalty}
\]

必须保持：

\[
\text{molecular support}
\ne
\text{enzyme abundance}
\ne
\text{enzyme capacity}
\ne
\text{flux}
\]

在当前小鼠肺癌脑转移与 SLC22A17 cKO 数据背景下，推荐将结果解释为：

> cell-type- and sample-aware multiome-supported metabolic hypotheses

并进一步通过培养基扰动、摄取/分泌实验、代谢组、蛋白组和同位素示踪验证。
