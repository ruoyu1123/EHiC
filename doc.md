# hicreate 算法说明

## 1. 软件设计目标

`hicreate` 是一个用于模�?Hi-C paired-end reads 的命令行程序。软件的目标是在给定参考基因组�?Hi-C 接触矩阵的条件下，生成符合输入接触频率分布的 150 bp paired-end FASTQ 数据。该程序既支持用户输入实验或人工构建�?Hi-C 矩阵，也支持在没有矩阵时根据参考基因组自动生成合成接触矩阵�?
软件的核心思想是：先将基因组划分为固定大小�?bin，并�?Hi-C 矩阵解释�?bin �?bin 之间的接触概率分布；随后通过限制性酶切模型建�?bin �?restriction fragment 之间的映射；最后按照矩阵权重采样虚拟连接事件，并从对应 restriction fragment 的末端生�?paired-end reads�?
因此，`hicreate` 的模拟过程不是简单地从全基因组均匀抽取 reads，而是�?Hi-C 矩阵为概率模型，显式模拟“接触矩�?-> 限制性片段连�?-> paired-end 测序 reads”的过程�?
## 2. 输入与输�?
软件主要输入包括�?
- 参考基因组 FASTA 文件�?- bin size，用于将基因组划分为固定大小�?genomic bins�?- 可选的 Hi-C 接触矩阵�?- 可选的 offset 文件，用于描述矩�?bin 与参�?contig/chromosome 的对应关系；
- 限制性内切酶识别序列�?- reads 数量、随机种子、线程数等模拟参数�?
软件输出为两�?FASTQ 文件�?
```text
<prefix>_R1.fastq
<prefix>_R2.fastq
```

每一�?reads 对应一个模拟的 Hi-C 连接事件。FASTQ header 中包含该 read pair 来源的两�?bin 编号，因此可以直接从输出 reads 重建模拟后的 Hi-C contact map。例如：

```text
@hicreate_123_bin10_bin45/1
```

表示�?123 �?read pair 来自 bin 10 �?bin 45 之间的接触事件�?
## 3. 总体算法流程

`hicreate` 的整体流程如下：

1. 读取参考基因组 FASTA�?2. 根据 bin size 构建全局 bin 坐标系统�?3. 读取用户输入�?Hi-C 矩阵，或生成合成 Hi-C 矩阵�?4. 根据 offset 将输入矩阵映射到当前参考基因组�?5. 对参考基因组进行 in silico restriction digestion�?6. 建立 bin �?restriction fragments 的索引；
7. 根据 contact matrix 权重采样 bin-bin 接触�?8. 从对�?bins 中采�?restriction fragments�?9. 构造虚�?ligation junction�?10. 生成 150 bp paired-end read templates�?11. 添加测序质量值和碱基替换错误�?12. 写出 R1/R2 FASTQ 文件�?
## 4. 参考基因组读取与分�?
程序首先读取输入 FASTA 文件。每条序列记录被解析为一�?contig，序列中的碱基统一转换为大写。空白字符会被忽略，FASTA header 的第一个字段作�?contig 名称�?
对于每条 contig，程序根据用户指定的 `bin_size` 计算�?bin 数：

```text
contig_bins = ceil(contig_length / bin_size)
```

多个 contig 会被拼接到同一个全局 bin 坐标系统中。程序为参考基因组生成 offset 表：

```text
contig    start_bin    end_bin
```

其中 `start_bin` 为该 contig 在全局 bin 坐标中的起始位置，`end_bin` 为开区间终点。该 offset 表用于判断两�?bins 是否属于同一 contig，也用于后续矩阵重映射和 cis/trans 接触统计�?
## 5. Hi-C 接触矩阵处理

### 5.1 自定义矩阵输�?
如果用户提供 `--matrix`，程序读取用户指定的 Hi-C 接触矩阵。矩阵支持两种格式�?
sparse 格式�?
```text
bin1    bin2    value
0       0       100
0       1       45
1       5       12
```

dense 格式�?
```text
100 45 12
45  90 18
12  18 80
```

矩阵格式可通过命令行参数显式指定：

```text
--matrix-format auto|sparse|dense
```

默认 `auto` 会根据输入内容自动判断。对�?3 x 3 dense 矩阵，由于其行数�?sparse 三列表达形式存在歧义，建议显式使用：

```text
--matrix-format dense
```

程序内部统一将矩阵转换为 sparse contact 列表�?
```text
Contact = (bin1, bin2, weight)
```

其中 `weight` 表示两个 bins 之间的相对接触强度。非正权重和非有限数值会被忽略。对�?sparse 输入，如果同一 bin pair 出现多次，程序会将其权重累加�?
### 5.2 输入矩阵�?cis/trans 含义

对于自定义输入矩阵，程序默认保留矩阵本身�?cis/trans 权重分布。也就是说，如果输入矩阵�?trans 接触较强，模�?reads �?trans 接触也会相应增强；如果输入矩阵中没有 trans 接触，程序不会默认强行添�?trans 接触�?
只有当用户显式指定：

```text
--trans-ratio X
```

时，程序才会对矩阵进�?cis/trans 重平衡，�?trans 权重占总权重的比例接近 `X`�?
```text
trans_fraction = trans_weight / (cis_weight + trans_weight)
```

该设计使得用户输入的真实或外部模拟矩阵能够作为主要概率模型，同时仍允许用户在需要时手动控制 trans 接触比例�?
### 5.3 矩阵重映�?
输入矩阵�?bin 坐标可能来自与当前参考基因组不同�?assembly、不同长度的模板，或不同�?contig 集合。因此，程序支持通过 offset 文件将输入矩阵映射到目标参考基因组�?
�?source offset �?target reference offset 中存在同�?contig 时，程序�?contig 内相对位置进行缩放映射：

```text
source_bin within source_contig
    -> proportional target_bin within target_contig
```

如果某些 contact 无法通过 contig 名称匹配，程序退化为全局比例缩放，即根据 source matrix 的�?bin 数和 target reference 的�?bin 数进行映射�?
如果输入矩阵只覆盖目标参考的一部分，程序会根据覆盖比例混合输入矩阵和合成背景矩阵。覆盖比例越高，输入矩阵权重占比越大；覆盖比例不足时，合成矩阵用于填补缺失信号�?
### 5.4 经验矩阵模型�?
仅依靠程序化规则从零生成 Hi-C 矩阵时，得到的热图形态可能难以稳定接近真实数据。为此，`hicreate` 支持经验矩阵模型库。模型库中的每个模型由真实或预先整理好的 Hi-C matrix �?offset 组成�?
```text
models/<model-name>/
  manifest.tsv
  matrix.tsv
  offset.tsv
```

`manifest.tsv` 记录矩阵文件、offset 文件和格式：

```text
name    human_cell_40mb
matrix  matrix.tsv
offset  offset.tsv
format  dense
```

支持的矩阵格式包�?text sparse、text dense 和紧�?binary sparse。用户可以通过�?
```text
--empirical-model human_cell_40mb
```

直接调用模型库�?
如果用户没有提供自己的矩阵，经验模型会被重映射到目标参考基因组，并作为完整 contact matrix 使用。如果用户同时提供了自己的矩阵和经验模型，则程序执行 trans-only replacement�?
1. 用户输入矩阵经过 remapping 后作�?base matrix�?2. base matrix 中所�?cis contacts 保留不变�?3. 经验模型矩阵经过 remapping 后只提取 trans contacts�?4. 经验模型 trans contacts resize、归一化后替换 base matrix 中原�?trans contacts�?5. 如果显式指定 `--trans-ratio`，替换后�?trans 总量按目标比例缩放；否则默认缩放到用户输入矩阵原�?trans 总量�?
这样可以实现�?
```text
final_matrix = user_cis + empirical_model_trans
```

该模式适合“用户热图提供样本特�?cis 结构，模型库提供可信 trans 架构”的场景，比完全手写 rule-based trans 模型更接近真�?Hi-C 数据�?
## 6. 合成 Hi-C 矩阵生成

当用户未提供输入矩阵时，`hicreate` 会根据参考基因组自动生成合成 Hi-C 矩阵。合成矩阵由 cis �?trans 两部分组成�?
### 6.1 cis 接触模型

cis 接触表示同一 contig/chromosome 内部的接触。程序使用距离衰减模型描�?cis 接触概率�?
```text
weight(d) �?1 / (d + 1)^alpha
```

其中 `d` 为两�?bins 之间的距离，`alpha` �?cis decay exponent。距离越近，接触权重越高。主对角线附近的 bins 因距离较短而具有较强接触信号�?
程序还支持设置最小和最�?cis 距离�?
```text
```

默认情况下允�?same-bin contact，从而保�?Hi-C 矩阵中常见的强主对角线�?
### 6.2 trans 接触模型

trans 接触表示不同 contigs/chromosomes 之间的接触。程序通过 `trans_ratio` 控制 trans 接触总权重占比，并结合物种模型和染色体空间排布模型生�?trans 接触�?
支持的空间模型包括：

- territory：染色体领地模型�?- rabl：Rabl 构型�?- rosette：Rosette-like 构型�?- nonrabl：非 Rabl 构型�?- custom：由用户输入的染色体三维空间布局�?
支持�?trans interaction 模型包括�?
- random：随机碰撞背景；
- territory / spatial / global-folding：基于染色体空间距离的接触；
- telomere：端�?亚端粒富集；
- centromere：着丝粒或近着丝粒区域富集�?- compartment / checkerboard：跨染色�?compartment-like 接触�?- hubs：若�?trans hotspot�?
这些模型用于在没有实验矩阵时生成具有一定生物学结构的合�?Hi-C contact map�?
### 6.3 用户指定染色体空间排�?
布局文件格式为：

```text
contig    x      y      z
chr1      0.0    0.0    0.0
chr2      0.8    0.1    0.2
chr3     -0.4    0.7   -0.1
```

其中 `contig` 必须�?FASTA �?offset 中的 contig 名称一致，`x/y/z` 为任意三维坐标系中的染色体中心位置。坐标单位不要求具有绝对物理意义，但相对距离会影�?trans 接触强度�?
当提供布局文件时，程序首先仍根据物种或 arrangement preset 生成一组默认染色体中心；随后用布局文件中匹配到�?contig 坐标覆盖默认中心。对于每一对不同染色体 `i` �?`j`，程序计算其三维距离�?
```text
d(i, j) = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
```

并将空间距离转换�?trans 染色体对权重�?
```text
arranged_weight(i, j) = 1 / (d(i, j) + 0.25)^gamma
```

同时，程序保留一个随机碰撞背景项�?
```text
```

最�?trans 染色体对的采样权重为�?
```text
pair_weight(i, j) = bins_i * bins_j * mixed_weight
```

其中 `bins_i` �?`bins_j` 分别表示两条染色体的 bin 数。这样，空间距离更近的染色体对会产生更强�?trans 接触信号；距离更远的染色体对则被削弱。该设计使用户能够显式模�?Rabl-like、rosette-like、territory-like、全局折叠或由外部 3D 建模得到的染色体空间构型�?

## 7. 限制性酶切模�?
Hi-C 实验中，基因�?DNA 通常经过限制性内切酶切割。因此，`hicreate` �?reads 生成前对参考基因组执行 in silico restriction digestion�?
用户可通过以下参数指定酶切识别序列�?
```text
--enzyme-site AAGCTT
```

程序会将输入的酶切位点统一转换为大写，因此 `aagctt` �?`AAGCTT` 等价�?
对于常见酶切位点，程序使用内�?cut offset�?
```text
AAGCTT  HindIII, A|AGCTT
GATC    DpnII/MboI, ^GATC
```

对于其他酶切位点，程序默认使�?motif 中点作为切割位置�?
每个 restriction fragment 记录如下信息�?
```text
fragment_id
contig_index
start
end
```

随后程序建立 bin �?fragments 的索引。如果一�?fragment 跨越多个 bins，则�?fragment 会被加入所有对�?bins 的候�?fragment 列表�?
## 8. 接触事件采样

程序�?contact matrix 中的每个 contact 权重解释为采样概率。对�?contact `i`�?
```text
P(contact_i) = weight_i / sum(weight)
```

程序首先根据所�?contacts 的权重构建累计分布，然后对每�?read pair 采样一�?bin-bin 接触�?
```text
(bin1, bin2)
```

接着，程序分别从 `bin1` �?`bin2` 对应�?restriction fragment 候选集合中随机选择一�?fragment，形成一个虚�?ligation event�?
如果 contact 来自同一�?bin 或相�?bins，两个端点可能采样到同一�?fragment。程序会尝试避免完全相同�?fragment 自连接，但当候�?fragment 数量有限时，仍允许这种事件存在。这种设计保留了�?bin、低酶切密度或局部强接触区域中的真实边界情况�?
## 9. 虚拟连接分子构建

对于每个采样到的 ligation event，程序从两个 restriction fragments 的末端向外截取序列，并构�?read template�?
每个 fragment 的末端方向会随机选择�?
- left end�?- right end 的反向互补序列�?
程序根据酶切位点构�?ligation junction。随后生�?read1 �?read2 的模板：

```text
read1_template = left_fragment_end + junction + right_fragment_end
read2_template = right_fragment_end + junction + left_fragment_end
```

如果 fragment 末端序列不足以生�?150 bp reads，程序会继续使用 junction 和另一�?fragment 序列补足。如果仍不足，则使用 `N` padding�?
该方法避免了显式构建完整的长 ligation molecule，同时保留了 Hi-C reads 由限制性片段末端和连接位点组成的核心特征�?
## 10. 测序质量与错误模�?
`hicreate` 在生�?read template 后，会进一步模�?Illumina-like 测序质量和碱基替换错误�?
质量值沿 read 位置变化：read 5' 端质量较高，�?3' 端逐渐下降。对于每个碱基位置，程序采样一�?Phred quality score `Q`，并根据 Phred 定义计算错误概率�?
```text
P(error) = 10^(-Q / 10)
```

当某个碱基发生错误时，程序会从其他三种碱基中随机选择一个替换。`N` 碱基不会被替换�?
最终，程序为每�?read 输出�?
```text
@read_name
sequence
+
quality
```

并分别写�?R1 �?R2 FASTQ 文件�?
## 11. 多线程与可复现�?
程序支持多线程生�?reads。为避免多线程写文件带来的锁竞争，`hicreate` �?reads 生成划分为多个固定大小的 block。worker threads 负责生成 FASTQ block，主线程�?block 顺序写入文件�?
这种设计具有两个优点�?
1. 内存占用不会随�?reads 数线性增长；
2. 输出顺序稳定，便于复现和下游比较�?
随机数由用户指定�?seed 初始化。每�?block 使用由全局 seed 派生的确定性随机流，从而减少线程调度对结果的影响�?
## 12. �?reads 重建 Hi-C 热图

由于每个 FASTQ header 中包含来�?bin 信息，输�?reads 可直接用于重建模�?contact map。对于每�?R1 header�?
```text
@hicreate_<id>_bin<i>_bin<j>/1
```

解析其中�?`i` �?`j`，并将矩阵中对应位置加一�?
```text
M[i, j] += 1
M[j, i] += 1, if i != j
```

这样得到�?reads-derived contact map 可以与输入矩阵或目标热图进行比较，用于评估模�?reads 是否保留了原�?Hi-C 矩阵结构�?
## 13. 算法特点

`hicreate` 的主要特点包括：

- �?Hi-C contact matrix 作为核心概率模型�?- 支持自定�?sparse/dense 矩阵输入�?- 支持 source matrix �?target reference �?bin-level 重映射；
- 保留 `--species-model` 作为物种级默认预设，同时允许用户用空间布局覆盖预设�?- 显式模拟限制性酶切片段；
- 根据 restriction fragment 端点构造虚�?ligation reads�?- 支持 Illumina-like 质量值和替换错误�?- 支持多线程流�?FASTQ 输出�?- 输出 read header 保留 contact bin 信息，便于质量控制和热图重建�?- 在自定义矩阵模式下默认保留输入矩阵的 cis/trans 权重比例，仅在用户显式指�?`--trans-ratio` 时进行重平衡�?
## 14. 方法学总结

总体而言，`hicreate` �?Hi-C reads 模拟问题分解为三个层次：

1. **接触概率�?*：由输入矩阵或合成矩阵定�?bin-bin 接触概率�?2. **基因组片段层**：由参考基因组和酶切位点定�?restriction fragments�?3. **测序 reads �?*：由 sampled ligation events 生成 paired-end reads，并加入测序质量和错误�?
这种分层设计使得软件既可以忠实地根据用户提供�?Hi-C 矩阵生成 reads，也可以在没有实验矩阵时生成具有基本生物学结构的合成 Hi-C reads。新增的染色体空间布局参数使模拟不再只能依赖粗粒度物种标签，而是可以直接围绕 Rabl、rosette、territory、checkerboard 或外�?3D 建模结果等具体空间构型设计实验。由于矩阵、offset、空间布局、酶切模型和随机种子均可由用户控制，`hicreate` 适用�?Hi-C pipeline 测试、算�?benchmark、参数敏感性分析以及不同接触模式下的模拟实验�?

## ��ǰ�������ģ��˵��

��ǰ�汾���Ƴ� rule-based Hi-C �������ɲ����Ͷ�Ӧ����������Ĭ������£�����û����ṩ --matrix���������� --species-model ѡ�� models/ �еľ������ģ�ͣ������� remap ��Ŀ��ο����������Ϊ���� contact matrix ʹ�á����û��ṩ�Զ�����������Զ������� cis �Ӵ��ṹ�����þ���ģ���е� trans �Ӵ����� resize��remap �͹�һ�����滻ԭʼ trans ���顣--trans-ratio �������������о��� trans �Ӵ�����������������ȱʧ�� rule-based trans contacts��


