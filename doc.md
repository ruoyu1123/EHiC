# hicreate 方法与算法说明

## 1. 软件目标

`hicreate` 根据参考基因组、限制性内切酶位点和接触矩阵模拟染色质构象捕获测序数据。当前支持三种 assay：

- `hic`：150 bp Illumina 风格双端 Hi-C reads。
- `porec`：Oxford Nanopore 风格的单分子多片段 Pore-C reads。
- `cifi`：PacBio HiFi 风格的单分子多片段 CiFi reads。

三种模式共用参考基因组、matrix、offset 和染色体级 matrix remap 流程，但 reads 构建与测序错误分别由独立模块实现。Hi-C 代码位于 `fragmenter.cpp` 和 `simulator.cpp`；Pore-C、CiFi 分别由 `porec.cpp`、`cifi.cpp` 调用 `long_read_common.cpp` 中的长 concatemer 生成框架。

基本命令格式为：

```bash
./hicreate ref.fa BIN_SIZE --assay hic|porec|cifi [options]
```

不指定 `--assay` 时使用 `hic`，因此旧的 Hi-C 命令保持兼容。

## 2. 共同输入和矩阵处理

### 2.1 参考基因组与全局 bin

程序读取 FASTA 中的每条 contig，以 header 第一个字段作为 contig 名称，并根据 `BIN_SIZE` 计算每条 contig 的 bin 数：

```text
contig_bins = ceil(contig_length / BIN_SIZE)
```

不同 contig 的 bin 按 FASTA 顺序连接成全局 bin 坐标。内部 reference offset 使用以下形式：

```text
contig  start_bin  end_bin
```

其中 `start_bin` 为闭区间起点，`end_bin` 为开区间终点。

### 2.2 Matrix 和 offset

Matrix 支持 sparse、dense 和 binary sparse 格式。Sparse 格式为：

```text
bin1  bin2  weight
```

每个正的有限 `weight` 被解释为对应 bin pair 的相对接触强度。用于直接采样时：

```text
P(contact_i) = weight_i / sum(all usable weights)
```

用户 offset 支持两种格式：

```text
contig  start_bin  end_bin
```

或 hicmap 四列格式：

```text
name  bin_offset  length  bin_num
```

四列格式按以下方式转换：

```text
start_bin = bin_offset
end_bin   = bin_offset + bin_num
```

`length` 列仅用于兼容，不参与 bin 区间计算。

### 2.3 Matrix 与参考基因组不匹配

程序首先检查 source offset 是否与 reference offset 完全一致。只有 contig 数量、名称、顺序和 bin 区间全部相同时，matrix 才直接进入模拟。

如果不一致，则执行染色体级 remap：

1. 优先匹配同名 source 和 target contig。
2. 剩余 contig 按 offset 顺序一一对应。
3. 每条 source contig 的完整 bin 区间按比例缩放到对应 target contig，而不是截断或用其他染色体补齐长度。
4. Cis contact 在同一对 source/target contig 内缩放。
5. Trans contact 按对应的 source contig pair 映射到 target contig pair。
6. 一个 source bin 映射到非整数 target 位置时，权重在相邻 target bins 之间分配。

默认要求 source offset 的 contig 数量不少于参考基因组。使用 `--force-contig-reuse` 可以强制复用 source contig；复制 target 之间的 trans block 从真实 source trans block 合成，并对最终 trans 总量归一化。

未提供 matrix 时，程序加载 `--species-model` 或 `--empirical-model` 指定的经验矩阵，并执行相同的 offset 检查和 remap。

同时提供用户 matrix 和 `--empirical-model` 时：

```text
final_matrix = user_matrix_cis + empirical_model_trans
```

经验 trans 会先映射到目标参考，再缩放到用户 matrix 原有的 trans 总量；显式指定 `--trans-ratio` 时，改为缩放到目标 trans 比例。

## 3. Hi-C reads 模拟方法

### 3.1 使用方式和输出

```bash
./hicreate ref.fa 1000 --assay hic \
  --matrix matrix.tsv --offset offset.tsv \
  --pairs 100000 --enzyme-site AAGCTT \
  --output-prefix sim
```

输出为：

```text
sim_R1.fastq
sim_R2.fastq
```

Hi-C 固定生成 150 bp paired-end reads。`--pairs N` 精确控制 R1、R2 各自的 FASTQ record 数。使用 coverage 时：

```text
pairs = ceil(coverage * reference_bases / (2 * 150))
```

### 3.2 In silico 酶切

Hi-C 默认酶切 motif 为 HindIII `AAGCTT`，切点为 `A|AGCTT`。`GATC` 按 DpnII/MboI 的 `^GATC` 处理；其他 motif 使用中点作为近似切点。

程序扫描每条 contig 并建立 restriction fragments：

```text
fragment_id, contig_index, start, end
```

如果一个 fragment 跨越多个 genomic bins，它会进入每个重叠 bin 的候选 fragment 集合。

### 3.3 接触和片段抽样

每个 read pair 首先按 matrix 权重抽取一个 `(bin1, bin2)` contact。随后分别从两个 bin 的候选 restriction fragments 中随机抽取一个片段。

程序会尝试避免两个端点选择同一个 fragment；如果候选片段不足，仍可能保留自连接事件。两个 fragment 分别随机选择左端或右端，右端以反向互补序列表示。

### 3.4 虚拟连接和随机 read 起点

程序为两个 fragment ends 构造两个方向的虚拟连接 scaffold：

```text
R1 scaffold = left_end  + ligation_junction + right_end
R2 scaffold = right_end + ligation_junction + left_end
```

随后在 scaffold 可用范围内随机抽取一个起始 offset，并从 R1、R2 scaffold 的对应位置分别截取 150 bp。因此 reads 并非总是固定从 fragment 的同一碱基开始。

如果 scaffold 不足 150 bp，剩余位置使用 `N` 补齐。

### 3.5 Illumina 质量和错误

质量值使用位置相关的 Illumina-like 模型：5' 端平均质量较高，向 3' 端逐渐下降。每个位置根据 Phred Q 计算 substitution 概率：

```text
P(error) = 10^(-Q / 10)
```

发生错误时，从其余三种碱基中随机选择替换碱基。当前 Hi-C 模块不模拟 indel。

## 4. Pore-C reads 模拟方法

### 4.1 建库过程对应关系

实验 Pore-C 的核心步骤是交联染色质、限制性酶切、近距离连接、解除交联、长片段 size selection、Nanopore adapter ligation 和单分子测序。一个 basecalled read 可以包含两个或更多来自不同基因组位置的 restriction fragments，因此能够表示 multi-way contact。

本程序模拟的是 Pore-C 的 basecalled FASTQ 层，不模拟 Nanopore 电流信号、basecaller、DNA 修饰或甲基化信号。

### 4.2 使用方式和输出

```bash
./hicreate ref.fa 1000 --assay porec \
  --matrix matrix.tsv --offset offset.tsv \
  --reads 10000 --threads 8 \
  --output-prefix porec_sim
```

Pore-C 默认使用 DpnII/MboI motif `GATC`。输出为：

```text
porec_sim_porec.fastq
porec_sim_porec_truth.tsv
```

Truth TSV 每行对应一个 concatemer segment：

```text
read_name
segment_order
contig
start
end
strand
matrix_bin
template_start
template_end
```

基因组坐标使用 0-based、end-exclusive 形式；`template_start/end` 是加入测序错误之前的 concatemer template 坐标。

### 4.3 Pore-C concatemer 长度

每条 read 的目标 template 长度来自 lognormal 分布。默认参数为：

```text
mean length     = 10000 bp
minimum length = 3000 bp
maximum length = 100000 bp
log sigma       = 0.80
```

抽样结果被限制在 minimum 和 maximum 之间。可以通过以下参数覆盖：

```text
--long-read-mean
--long-read-min
--long-read-max
--long-read-sigma
```

`--long-read-sigma 0` 表示使用固定目标长度。

### 4.4 首个 contact 和 anchor 模型

输入 matrix 是二元接触矩阵，不能唯一决定三个及以上片段的联合分布。因此程序使用 anchor 模型构造高阶 concatemer。

第一步严格按 matrix 总权重抽取一个 contact：

```text
(bin1, bin2) ~ matrix weights
```

两个端点的顺序随机交换，交换后的第一个 bin 成为该 concatemer 的 `anchor_bin`。前两个 segment 分别来自这两个 bins。

对于第三个及后续 segment，程序从 anchor 的邻接接触分布中抽样：

```text
P(next_bin = j | anchor = i)
    = weight(i, j) / sum_k weight(i, k)
```

这里始终使用同一个 anchor，而不是从上一 segment 继续随机游走。这样可以避免一次 trans 跳转让后续大量片段持续停留在另一条染色体，进而在 all-vs-all 展开时人为放大 trans 信号。

### 4.5 Segment 和连接序列构造

对每个 sampled bin，程序从与该 bin 重叠的 restriction fragments 中随机选一个片段，并尽量避免连续两次使用同一 fragment。

每个 segment 的方向独立随机：

- `+`：使用参考序列正向片段。
- `-`：使用对应片段的反向互补序列。

相邻 segments 之间插入根据酶切 motif 构造的 canonical ligation junction。程序持续增加 segments，直到达到抽样的目标 template 长度或 `--max-segments`。默认最大值为 256；如果先达到 segment 上限，程序会报告 warning 和 `Segment-limited reads` 数量。

首尾 segment 必要时可以只使用 restriction fragment 的一部分，以模拟 size selection 或测序分子边界。每条输出 read 至少包含两个 truth segments。

### 4.6 Nanopore 风格测序错误

Pore-C 默认平均 QV 为 13，位置 QV 从标准差为 2 的正态分布抽样，并限制在 Q2 到 Q41。每个 template base 的总错误概率为：

```text
P(error) = 10^(-Q / 10)
```

错误类型默认按以下比例分配：

```text
deletion      45%
substitution  25%
insertion     30%
```

Deletion 不输出该 template base；substitution 输出随机替换碱基；insertion 在当前碱基后增加随机碱基。最终 FASTQ sequence 和 quality 长度始终一致。`--long-read-qv` 可以覆盖平均 QV。

## 5. CiFi reads 模拟方法

### 5.1 建库过程对应关系

CiFi 将 3C concatemer 与 PacBio HiFi 测序结合。实验流程包括交联、DpnII 或 HindIII 酶切、近距离连接、解除交联、whole-genome PCR enrichment、SMRTbell library preparation、长片段 size selection 和 CCS/HiFi sequencing。

本程序直接模拟最终 HiFi consensus FASTQ，不模拟 polymerase raw read、SMRTbell 环形序列、不同 pass 数量或 CCS consensus 计算过程。

### 5.2 使用方式和输出

```bash
./hicreate ref.fa 1000 --assay cifi \
  --matrix matrix.tsv --offset offset.tsv \
  --coverage 5 --threads 8 \
  --output-prefix cifi_sim
```

CiFi 默认使用 `GATC`，也可以显式指定 HindIII：

```text
--enzyme-site AAGCTT
```

输出为：

```text
cifi_sim_cifi.fastq
cifi_sim_cifi_truth.tsv
```

Truth TSV 的字段和坐标规则与 Pore-C 相同。

### 5.3 CiFi concatemer 构造

CiFi 与 Pore-C 使用相同的 restriction fragment index、首个 matrix contact 抽样和 anchor 条件分布。区别主要来自 read 长度、测序精度和 PCR duplicate 模型，而不是使用另一套 contact matrix。

默认长度参数来自已发表 CiFi 数据的近似统计：

```text
mean length     = 9350 bp
minimum length = 5000 bp
maximum length = 25000 bp
log sigma       = 0.65
```

限制性酶会自然影响 segments 数量。DpnII 位点更密集，同样长度的 concatemer 通常包含更多、较短的 segments；HindIII 位点更稀疏，通常包含更少、较长的 segments。程序不直接强制使用“DpnII 17 段”或“HindIII 2 段”，而是由参考序列中的真实 motif 密度和目标 read 长度共同决定。

### 5.4 HiFi 质量和错误

CiFi 默认平均 QV 为 38，QV 标准差为 1.5，并限制在 Q2 到 Q41。错误概率仍按 Phred 定义计算，但显著低于 Pore-C。

默认错误类型比例为：

```text
substitution  55%
deletion      25%
insertion     20%
```

这些错误施加在 concatemer template 上，输出代表完成 CCS 后的高准确度 HiFi read。

### 5.5 PCR template duplicate

CiFi 建库包含 PCR enrichment。程序默认以 1.8% 概率复用同一 FASTQ block 中先前生成的 concatemer template：

```text
template_duplicate_rate = 0.018
```

Duplicate read 使用相同的 source segments、坐标、方向和连接结构，但重新独立采样 HiFi 测序错误。Duplicate read 名称包含：

```text
_dup<SOURCE_READ_ID>
```

Pore-C 不使用该 PCR duplicate 模型。

## 6. 三种模式比较

| 特征 | Hi-C | Pore-C | CiFi |
|---|---|---|---|
| 输出形态 | PE150 | 单端长 read | 单端 HiFi 长 read |
| 每条 read 的 segment 数 | 2 个连接端 | 至少 2 个 | 至少 2 个 |
| 默认酶切 motif | `AAGCTT` | `GATC` | `GATC` |
| Matrix 用法 | 每个 pair 直接抽一个 contact | 首对直接抽样，后续围绕 anchor | 首对直接抽样，后续围绕 anchor |
| 默认平均长度 | 2 x 150 bp | 10 kb | 9.35 kb |
| 默认平均 QV | Illumina 位置模型 | Q13 | Q38 |
| Indel | 不模拟 | 模拟 | 模拟，比例较低 |
| PCR duplicate | 不模拟 | 不模拟 | 默认 1.8% template duplicate |
| Truth TSV | 无 | 有 | 有 |

## 7. Reads 数量和 coverage

Hi-C 使用：

```text
--pairs N
```

Pore-C、CiFi 使用：

```text
--reads N
```

`--reads` 精确控制单端长 read record 数量。长 read coverage 模式先按配置的平均长度换算 read 数：

```text
reads = ceil(coverage * reference_bases / configured_mean_read_length)
```

由于实际长度来自随机分布，并且 insertion/deletion 会改变输出碱基数，最终 coverage 可能略高或略低。程序结束时按实际 FASTQ 碱基数报告：

```text
actual_coverage = sequenced_bases / reference_bases
```

`--coverage` 不能与 `--pairs` 或 `--reads` 同时使用。

## 8. 多线程、内存和可复现性

三种模式都将 reads 分成有界 block，由 worker threads 负责抽样、构建模板、模拟错误和格式化 FASTQ，主线程按 block 编号顺序写盘。

长 read block 大小会根据平均 read 长度动态调整，避免超长 read 时一次预留过多内存。ready queue 有容量限制，但 writer 正在等待的 block 可以越过容量限制进入队列，从而避免乱序 worker 导致死锁。

每个 block 使用由全局 `--seed` 和 block 编号派生的随机流。因此在参数、输入和 seed 相同的情况下，改变线程数不会改变 FASTQ 和 truth TSV 内容。

## 9. 高阶接触的解释和限制

Pore-C 和 CiFi 的一条 concatemer 可以拆分成 `n` 个 segments，并产生 `n * (n - 1) / 2` 个 all-vs-all pairwise contacts。输入 matrix 只提供二元边缘分布，不包含以下信息：

- 哪些三个或更多 loci 同时存在于一个细胞复合物中。
- Multi-way complex 中各 segments 的联合概率。
- Segment 在连接分子中的真实物理顺序。

因此不存在从普通 Hi-C matrix 到唯一 Pore-C/CiFi concatemer 集合的精确逆变换。当前 anchor 模型能保持首个 contact 精确服从 matrix，并避免明显的跨染色体随机游走，但将所有 concatemer 做 all-vs-all 展开后，不保证每个像素或 cis/trans 比例与输入 matrix 完全相等。

需要严格评估时，应使用 truth TSV 重建：

1. 相邻 segment contacts。
2. Anchor-to-all contacts。
3. All-vs-all contacts。

比较三种统计与输入 matrix 后，再针对具体下游 pipeline 校准长度分布、`--max-segments`、酶切 motif 或额外的高阶模型。

## 10. 方法依据

Pore-C 模型参考：

- Deshpande et al. *Nanopore sequencing of DNA concatemers reveals higher-order features of chromatin structure*. DOI: `10.1101/833590`。
- Zhong et al. *High-throughput Pore-C reveals the single-allele topology and cell type-specificity of 3D genome folding*. DOI: `10.1038/s41467-023-36899-x`。

CiFi 模型参考：

- McGinty et al. *CiFi: accurate long-read chromosome conformation capture with low-input requirements*. DOI: `10.1038/s41467-025-66918-y`。
- *CiFi: 3C Library Preparation for PacBio HiFi Sequencing*. DOI: `10.17504/protocols.io.4r3l21zxpg1y/v1`。
