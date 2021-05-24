# set workpath and load R packages
setwd("/ifs1/Grp8/liuzhe/FengLab/gliomaFP/scRNA-seq/")
rm(list = ls())
#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
#remotes::install_github("jlmelville/uwot")
#install.packages('Seurat')
library("Seurat")
#install.packages('tidyverse')
library("tidyverse")
library("Matrix")
library("scales")
library("cowplot")
library("RCurl")
options(stringsAsFactors = F)
library("scRNAseq")
library("scater")
# devtools::install_github(repo = "satijalab/seurat", ref = "loom")
library("loomR")
library("patchwork")
library("scuttle")
library("dplyr")
library("tibble")
library("HCAData")
library("SingleR")
library("org.Hs.eg.db")
library("clusterProfiler")
library("vroom")
library("celldex")
library("dittoSeq")
set.seed(2020)
library("monocle")
suppressPackageStartupMessages(library("SingleCellExperiment"))
library("ggplot2"); theme_set(theme_bw())
library("DuoClustering2018")
require(scry)

# load the scRNA-seq data
# create each individual Seurat object for every sample
# create Seurat object gliomaF2, gliomaF3, and gliomaP1
for(file in c("gliomaP1","gliomaP2","gliomaP3","gliomaF2","gliomaF3","gliomaF4")){
	seurat_data<-Read10X(data.dir=paste0("1_counts/",file))
	seurat_obj<-CreateSeuratObject(counts=seurat_data,min.cells = 3,min.features=200,project=file)
	assign(file,seurat_obj)
}
gliomaP1<-NormalizeData(gliomaP1)
gliomaP2<-NormalizeData(gliomaP2)
gliomaP3<-NormalizeData(gliomaP3)
gliomaF2<-NormalizeData(gliomaF2)
gliomaF3<-NormalizeData(gliomaF3)
gliomaF4<-NormalizeData(gliomaF4)
gliomaP.normalized.combined <- merge(gliomaP1, y = c(gliomaP2, gliomaP3), add.cell.ids = c("P1", "P2", "P3"), project = "gliomaP", merge.data = TRUE)
gliomaF.normalized.combined <- merge(gliomaF2, y = c(gliomaF3, gliomaF4), add.cell.ids = c("F2", "F3", "F4"), project = "gliomaF", merge.data = TRUE)
GetAssayData(gliomaP.normalized.combined)[1:10, 1:15]
gliomaP.normalized.combined
GetAssayData(gliomaF.normalized.combined)[1:10, 1:15]
gliomaF.normalized.combined
pdf("2_scRNA-seq/figures/workflow/gliomaP.counts.vs.features.pdf")
plot(x=gliomaP.normalized.combined@meta.data$nCount_RNA,y=gliomaP.normalized.combined@meta.data$nFeature_RNA)
dev.off()
pdf("2_scRNA-seq/figures/workflow/gliomaF.counts.vs.features.pdf")
plot(x=gliomaF.normalized.combined@meta.data$nCount_RNA,y=gliomaF.normalized.combined@meta.data$nFeature_RNA)
dev.off()
# check the metadata in the new Seurat objects
head(gliomaP.normalized.combined@meta.data)
tail(gliomaP.normalized.combined@meta.data)
head(gliomaF.normalized.combined@meta.data)
tail(gliomaF.normalized.combined@meta.data)
# Create .RData object to load at any time
save(gliomaP.normalized.combined, file="2_scRNA-seq/data/gliomaP.combined.RData")
save(gliomaF.normalized.combined, file="2_scRNA-seq/data/gliomaF.combined.RData")
# Create a merged Seurat object
gliomaP.normalized.combined<-NormalizeData(gliomaP.normalized.combined)
gliomaF.normalized.combined<-NormalizeData(gliomaF.normalized.combined)
merged_seurat<-merge(gliomaP.normalized.combined,y=gliomaF.normalized.combined, project = "Glioma", add.cell.id=c("gliomaP","gliomaF"), merge.data = TRUE)
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
# Seurat会自动为每个细胞创建一些元数据：
#View(merged_seurat@meta.data)
# 添加的列
# orig.ident ：通常包含样本标识（如果已知），通常默认project为我们为其分配的身份
# nCount_RNA ：每个细胞的UMI数量
#nFeature_RNA ：每个细胞检测到的基因数量
# 我们需要计算一些用于绘图的其他指标：
# number of genes detected per UMI: 这个度量让我们对数据集的复杂性有了一个概念(每个UMI检测到的基因越多，我们的数据就越复杂)
# mitochondrial ratio: 这个度量将给我们提供一个源自线粒体基因的细胞读数的百分比。
# 每个细胞的每个UMI的基因数量非常容易计算，我们将对结果进行log10转换，以便更好地在样本之间进行比较。
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Seurat有一个方便的功能，可以让我们计算映射到线粒体基因的转录本的比例。PercentageFeatureSet()将采用一定模型并搜索基因标识符。此功能可以轻松计算属于每个细胞的可能功能子集的所有计数的百分比。这里的计算只是将属于该集合的要素的计数槽中存在的矩阵的列和除以所有要素的列和，然后乘以100。由于我们要绘制比率值，所以我们将反转这一步，然后除以100
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
# 注意：提供的模式（“ ^ MT-”）适用于人类基因名称。Mt适用于大鼠基因名称。 您可能需要根据您感兴趣的生物进行调整。如果您没有使用基因名称作为基因ID，则此功能将无法使用。我们有代码可用于自行计算该指标。
# 同时还需要将其他信息添加到QC指标的元数据中，例如单元ID、条件信息和各种指标。虽然使用$操作符将信息直接添加到Seurat对象的元数据槽非常容易，但是我们选择把数据框提取到一个单独的变量中。通过这种方式，我们可以继续插入QC分析所需的其他指标，而不会有影响merded_seurat对象的风险。
# 我们将通过从seurat对象提取meta.data槽来创建元数据数据框：
# Create metadata dataframe
metadata <- merged_seurat@meta.data
# 可以看到每个细胞ID都有一个gliomaP_P1，gliomaP_P2等等的前缀，正如我们在合并Seurat对象时指定的那样。这些前缀应该与Orig.ident中列出的样本匹配。让我们首先添加一个包含细胞ID的列，并更改当前列名以使其更直观：
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
# 现在，根据细胞前缀获取每个细胞的样本名称
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^gliomaP_"))] <- "gliomaP"
metadata$sample[which(str_detect(metadata$cells, "^gliomaF_"))] <- "gliomaF"
head(metadata)
tail(metadata)
# 现在，您已经设置好了评估数据质量所需的指标！最终的元数据表将包含对应于每个细胞的行，以及包含有关这些细胞的信息的列：
# 将更新的元数据保存到我们的Seurat对象
# 在评估指标之前，我们可以把迄今为止完成的所有工作保存回Seurat对象中，这样方便以后调用。我们只需将数据框分配到meta.data插槽即可完成此操作。
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
# Create .RData object to load at any time
save(merged_seurat, file="2_scRNA-seq/data/merged_filtered_seurat.RData")
# 评估质量指标
# 现在我们已经生成了要评估的各种指标，我们可以通过可视化来探索它们。我们将评估各种指标，然后决定哪些细胞质量较低，应该从分析中删除：
# 细胞计数
# 每个细胞的UMI计数
# 每个细胞检测到的基因
# UMI与检测到的基因
# 线粒体比率
# Novelty
# What about doublets?在单细胞RNA测序实验中，双胞体是由两个细胞产生的。它们通常是由于细胞分选或捕获中的错误引起的，特别是在涉及数千个细胞的基于液滴的协议中。当目标是描述单细胞水平的群体特征时，双峰显然是不可取的。具体地说，他们可能错误地暗示存在实际并不存在的中间群体或短暂状态。因此，需要移除双峰文库，以便它们不会影响结果的解释。
# Why aren’t we checking for doublets? 许多工作流程使用UMI或基因的最大阈值，其想法是检测到的大量读数或基因表明存在多个细胞。尽管此理由似乎很直观，但并不准确。同样，许多用于检测双峰的工具倾向于去除具有中间或连续表型的细胞，尽管它们在具有非常离散的细胞类型的数据集上可能会很好地工作。Scrublet是双态检测的流行工具，但我们尚未对其进行充分的基准测试。目前，我们建议您此时不设置任何阈值。当我们确定了每个簇的标记时，建议您探索这些标记，以确定这些标记是否适用于一种以上的细胞类型。
# Cell counts
# 细胞计数由检测到的唯一细胞条形码的数量确定。对于此实验，预计将有12,000 -13,000个细胞。
# 在理想的情况下，您希望唯一的细胞条形码的数量与您加载的细胞的数量相对应。但是，情况并非如此，因为细胞的捕获率仅是所装载细胞的比例的一部分。例如，与10x相比inDrops的细胞捕获效率更高(70-80%)，而10x的效率可以下降到50-60%。
# 注意：如果用于文库制备的细胞浓度不准确，捕获效率可能会显得很低。细胞浓度不应由FACS机器或生物分析仪测定(这些工具对浓度测定不准确)，而应使用血细胞计数仪或自动细胞计数器计算细胞浓度。
# 细胞编号也可能因方法而异，产生的细胞编号比我们加载的要高得多。例如，在inDrops方法中，细胞条形码存在于水凝胶中，并与单个细胞和裂解/反应混合物封装在液滴中。虽然每个水凝胶都应该有一个与之相关的细胞条形码，但有时一个水凝胶可以有多个细胞条形码。类似地，使用10X协议，有可能只获得乳液液滴(GEM)中的条形码珠子，而没有实际的细胞。这两种情况，加上死亡细胞的存在，都可能导致比细胞更多的细胞条形码。
# Visualize the number of cell counts per sample
pdf("2_scRNA-seq/figures/pre-operation/thenoofcellcountspersample.pdf")
metadata %>% 
      ggplot(aes(x=sample, fill=sample)) + 
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      ggtitle("NCells")
dev.off()
# 我们看到每个样本都超过15，000个细胞，这比预期的12-13,000个要多得多。很明显，我们现在可能有一些垃圾“细胞”。
# UMI counts (transcripts) per cell
# 每个细胞的UMI计数通常应高于500，这是我们预期的下限。如果UMI计数在500-1000计数之间，则可以使用，但可能应该对细胞进行更深的测序。
# Visualize the number UMIs/transcripts per cell
pdf("2_scRNA-seq/figures/pre-operation/thenoUMIspercell.pdf")
metadata %>% 
      ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      theme_classic() +
      ylab("Cell density") +
      geom_vline(xintercept = 500)
dev.off()
# 我们可以看到，两个样本中的大多数单元都具有1000个UMI或更高，这非常好。
# Genes detected per cell
# 我们对基因检测的期望值与UMI检测的期望值相似，尽管可能比UMIs略低。对于高质量数据，比例直方图应包含表示封装的细胞的单个大峰值。如果我们看到主峰右边的一个小肩膀(我们的数据中没有出现)，或者细胞的双峰分布，这可能表明了一些事情。可能是因为某些原因导致一组细胞发生错误。也可能是存在生物上不同类型的细胞(即，静止的细胞群体、不太复杂的目标细胞)，还有就是可能是一种类型比另一种小得多(即，高计数的细胞可能是尺寸较大的细胞)。因此，应该使用我们在本课中描述的其他指标来评估此阈值。
#UMIs vs. genes detected
#通常一起评估的两个指标是UMI的数量和每个细胞检测到的基因数量。在这里，我们绘制了基因数量与线粒体reads所占比例UMI数量的关系图。线粒体reads部分只在检测到很少基因的特别低计数的细胞中才高(浅蓝色)。这可能是损伤/死亡的细胞，其细胞质的mRNA已经通过破裂的膜泄漏出来，因此，只有位于线粒体的mRNA仍然是保守的。这些细胞被我们的计数和基因数量阈值过滤掉。联合可视化计数和基因阈值可显示联合过滤效果。
#质量差的细胞很可能每个细胞的基因和UMI都很低，并且与图左下象限的数据点相对应。好的细胞通常会表现为每个细胞有更多的基因和更高数量的UMI。
#通过此图，我们还评估了线的斜率，以及图的右下角象限中数据点的任何散布情况。这些细胞有大量的UMI，但只有几个基因。这些可能是濒临死亡的细胞，但也可能代表一个低复杂性细胞类型的群体(即红细胞)。
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("2_scRNA-seq/figures/pre-operation/thecorrelationbetweengenesdetectedandnoofUMIs.pdf")
metadata %>% 
      ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
      geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
      stat_smooth(method=lm) +
      scale_x_log10() + 
      scale_y_log10() + 
      theme_classic() +
      geom_vline(xintercept = 500) +
      geom_hline(yintercept = 250) +
      facet_wrap(~sample)
dev.off()
# Mitochondrial counts ratio
# 这一指标可以确定死亡或濒临死亡的细胞是否存在大量的线粒体污染。我们将线粒体计数质量差的样品定义为超过0.2线粒体比率标记的细胞，除非您希望样品中有这种情况。
# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("2_scRNA-seq/figures/pre-operation/theMitochondrialcountsratio.pdf")
metadata %>% 
      ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      theme_classic() +
      geom_vline(xintercept = 0.2)
dev.off()
# Complexity
# 我们可以看到，我们对每个细胞测序较少的样本具有更高的整体复杂性，这是因为我们还没有开始对这些样本的任何给定基因进行饱和测序。这些样本中的异常值细胞可能是RNA种类比其他细胞简单的细胞。有时，我们可以通过这个指标检测到低复杂性细胞类型的污染，比如红细胞。一般来说，我们预计novelty得分在0.80以上。
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
pdf("2_scRNA-seq/figures/pre-operation/theoverallcomplexityofthegeneexpression.pdf")
metadata %>%
      ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8)
dev.off()
# 注意：Reads per cell是另一个值得研究的指标；但是，使用的工作流需要保存此信息以供评估。通常，可以使用此度量标准来查看所有样本，每个样本的峰值在相对相同的位置，每个细胞的读数介于10,000和100,000之间。
# 过滤
# 总之，孤立地考虑这些质量控制指标中的任何一个都可能导致对细胞信号的误解。例如，线粒体计数相对较高的细胞可能参与呼吸过程，可能是您想要保留的细胞。同样，其他指标也可以有其他生物学解释。因此，在设置阈值时，请始终考虑这些指标的共同影响，并将其设置为尽可能宽松，以避免无意中过滤掉可行的细胞群体。
# Cell-level filtering
# 现在我们已经可视化了各种指标，我们可以决定要应用的阈值，这将导致删除低质量的细胞。通常，前面提到的建议是一个粗略的指导原则，而具体的实验需要设定所选的确切阈值。我们将使用以下阈值：
# nUMI > 500
# nGene > 250
# log10GenesPerUMI > 0.8
# mitoRatio < 0.2
#进行过滤，使用以下subset()函数：
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
# Gene-level filtering
# 在我们的数据中，我们将有许多零计数的基因。这些基因可以极大地降低细胞的平均表达量，所以我们将把它们从我们的数据中删除。首先，我们将删除所有细胞中零表达的基因。此外，我们还将根据prevalence执行一些过滤。如果一个基因仅在少数细胞中表达，那么它就没有什么特别的意义，因为它仍会降低未在其中表达的所有其他细胞的平均值。对于我们的数据，我们选择只保留在10个或更多细胞中表达的基因。
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# 重新评估QC指标
# 执行过滤之后，建议回顾一下指标，以确保您的数据符合您的期望并适合进行下游分析。
# 保存过滤后的细胞
# 基于这些QC指标，我们将识别出任何不合格的样本，并继续进行过滤后的单元格。通常，我们使用不同的过滤条件来迭代QC指标。它不一定是线性过程。满足过滤条件后，我们将保存过滤后的细胞对象以进行聚类和标记识别。
# Create .RData object to load at any time
save(filtered_seurat, file="2_scRNA-seq/figures/pre-operation/seurat_filtered.RData")
gliomaF<-gliomaF.normalized.combined
gliomaP<-gliomaP.normalized.combined
gliomaP$log10GenesPerUMI <- log10(gliomaP$nFeature_RNA) / log10(gliomaP$nCount_RNA)
gliomaP$mitoRatio <- PercentageFeatureSet(object = gliomaP, pattern = "^MT-")
gliomaP$mitoRatio <- gliomaP@meta.data$mitoRatio / 100
gliomaPmetadata <- gliomaP@meta.data
gliomaPmetadata$cells <- rownames(gliomaPmetadata)
gliomaPmetadata <- gliomaPmetadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
gliomaP
gliomaP@meta.data <- gliomaPmetadata
counts <- GetAssayData(object = gliomaP, slot = "counts")
gliomaP <- CreateSeuratObject(counts, meta.data = gliomaP@meta.data)
gliomaP$label <- "gliomaP"
gliomaP_norm <- NormalizeData(gliomaP, normalization.method = "LogNormalize", scale.factor = 10000)
gliomaP_norm <- FindVariableFeatures(gliomaP_norm, selection.method = "vst", nfeatures = 2000)
pdf("2_scRNA-seq/figures/pre-operation/gliomaP_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(gliomaP, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(gliomaP, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(gliomaP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("2_scRNA-seq/figures/pre-operation/gliomaP_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
top30 <- head(VariableFeatures(gliomaP_norm), 30)
pdf("2_scRNA-seq/figures/pre-operation/gliomaP_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(gliomaP_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
gliomaF$log10GenesPerUMI <- log10(gliomaF$nFeature_RNA) / log10(gliomaF$nCount_RNA)
gliomaF$mitoRatio <- PercentageFeatureSet(object = gliomaF, pattern = "^MT-")
gliomaF$mitoRatio <- gliomaF@meta.data$mitoRatio / 100
gliomaFmetadata <- gliomaF@meta.data
gliomaFmetadata$cells <- rownames(gliomaFmetadata)
gliomaFmetadata <- gliomaFmetadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
gliomaF
gliomaF@meta.data <- gliomaFmetadata
counts <- GetAssayData(object = gliomaF, slot = "counts")
gliomaF <- CreateSeuratObject(counts, meta.data = gliomaF@meta.data)
gliomaF$label <- "gliomaF"
gliomaF_norm <- NormalizeData(gliomaF, normalization.method = "LogNormalize", scale.factor = 10000)
gliomaF_norm <- FindVariableFeatures(gliomaF_norm, selection.method = "vst", nfeatures = 2000)
pdf("2_scRNA-seq/figures/pre-operation/gliomaF_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(gliomaF, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(gliomaF, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(gliomaF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("2_scRNA-seq/figures/pre-operation/gliomaF_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
top30 <- head(VariableFeatures(gliomaF_norm), 30)
pdf("2_scRNA-seq/figures/pre-operation/gliomaF_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(gliomaF_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
# filter gliomaP
filtered_gliomaP <- subset(x = gliomaP, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
counts <- GetAssayData(object = filtered_gliomaP, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_gliomaP <- CreateSeuratObject(filtered_counts, meta.data = gliomaP@meta.data)
filtered_gliomaP$label <- "gliomaP"
save(filtered_gliomaP, file="2_scRNA-seq/data/filtered_gliomaP.RData")
# filter gliomaF
filtered_gliomaF <- subset(x = gliomaF, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
counts <- GetAssayData(object = filtered_gliomaF, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_gliomaF <- CreateSeuratObject(filtered_counts, meta.data = filtered_gliomaF@meta.data)
filtered_gliomaF$label <- "gliomaF"
save(filtered_gliomaF, file="2_scRNA-seq/data/filtered_gliomaF.RData")
filtered_merged_seurat<-merge(filtered_gliomaP,y=filtered_gliomaF, project = "filtered_Glioma", add.cell.id=c("filtered_gliomaP","filtered_gliomaF"))
head(filtered_merged_seurat@meta.data)
tail(filtered_merged_seurat@meta.data)
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaF_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(filtered_gliomaF, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(filtered_gliomaF, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(filtered_gliomaF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaF_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaP_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(filtered_gliomaP, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(filtered_gliomaP, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(filtered_gliomaP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaP_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
filtered_gliomaF_norm <- NormalizeData(filtered_gliomaF, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_gliomaF_norm <- FindVariableFeatures(filtered_gliomaF_norm, selection.method = "vst", nfeatures = 2000)
save(filtered_gliomaF_norm, file="2_scRNA-seq/data/filtered_gliomaF_norm.RData")
top30 <- head(VariableFeatures(filtered_gliomaF_norm), 30)
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaF_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(filtered_gliomaF_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
filtered_gliomaP_norm <- NormalizeData(filtered_gliomaP, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_gliomaP_norm <- FindVariableFeatures(filtered_gliomaP_norm, selection.method = "vst", nfeatures = 2000)
save(filtered_gliomaP_norm, file="2_scRNA-seq/data/filtered_gliomaP_norm.RData")
top30 <- head(VariableFeatures(filtered_gliomaP_norm), 30)
pdf("2_scRNA-seq/figures/pre-operation/filtered_gliomaP_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(filtered_gliomaP_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
save(filtered_seurat, file="2_scRNA-seq/data/filtered_seurat20210510.RData")
glioma.anchors <- FindIntegrationAnchors(object.list = list(filtered_gliomaP_norm, filtered_gliomaF_norm), dims = 1:30)
save(glioma.anchors, file="2_scRNA-seq/data/integrated.anchors_seurat20210510.RData")
glioma.combined <- IntegrateData(anchorset = glioma.anchors, dims = 1:30)
glioma.combined <- FindVariableFeatures(glioma.combined, selection.method = "vst", nfeatures = 2000)
save(glioma.combined, file="2_scRNA-seq/data/integrated.combined_seurat20210510.RData")
DefaultAssay(glioma.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
glioma.combined <- ScaleData(glioma.combined, verbose = FALSE, vars.to.regress = c("nUMI", "mitoRatio"))
save(glioma.combined, file="2_scRNA-seq/data/glioma.combined_scaled.RData")
glioma.combined <- RunPCA(glioma.combined, npcs = 30, verbose = FALSE)
print(glioma.combined[["pca"]], dims = 1:5, nfeatures = 5)
pdf("2_scRNA-seq/figures/pre-operation/VizDimLoadings.pdf")
VizDimLoadings(glioma.combined, dims = 1:2, reduction = "pca")
dev.off()
pdf("2_scRNA-seq/figures/pre-operation/DimPlot.pdf")
DimPlot(glioma.combined, reduction = "pca")
dev.off()
pdf("2_scRNA-seq/figures/pre-operation/DimHeatmap.pc1.pdf")
DimHeatmap(glioma.combined, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("2_scRNA-seq/figures/pre-operation/DimHeatmap.all.pdf")
DimHeatmap(glioma.combined, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
glioma.combined <- JackStraw(glioma.combined, num.replicate = 100, dims = 30)
glioma.combined <- ScoreJackStraw(glioma.combined, dims = 1:30)
pdf("2_scRNA-seq/figures/pre-operation/Determine_dimensionality.pdf", width = 24, height = 18)
p1 <- JackStrawPlot(glioma.combined, dims = 1:30)
p2 <- ElbowPlot(glioma.combined,ndims = 30)
plot_grid(p1, p2)
dev.off()
save(glioma.combined, file="2_scRNA-seq/data/integrated.combined_beforepcs.RData")
# determine the resolution
#library(Seurat)
#library(clustree)
#sce <- FindClusters(
#  object = sce,
#  resolution = c(seq(.1,1.6,.2))
#)
#pdf("2_scRNA-seq/figures/pre-operation/clusters.pdf",width=30,height=15)
#clustree(sce@meta.data, prefix = "integrated_snn_res.")
#colnames(sce@meta.data)
#dev.off()

r=0.75
glioma.combined.pca29 <- FindNeighbors(glioma.combined, reduction = "pca", dims = 1:29)
glioma.combined.pca29 <- FindClusters(glioma.combined.pca29, resolution = r)
levels(glioma.combined.pca29)
save(glioma.combined.pca29,file="2_scRNA-seq/data/glioma.combined.pca29.res0.75.RData")
glioma.pca29.markers <- FindAllMarkers(object = glioma.combined.pca29, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
save(glioma.pca29.markers,file="2_scRNA-seq/data/glioma.pca29.markers.res0.75.RData")

###################################################################################################################################
# set workpath and load R packages
setwd("/home/liuzhe/FengLab/glioma_project/1_singlecellrnaseq")
rm(list = ls())
#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
#remotes::install_github("jlmelville/uwot")
#install.packages('Seurat')
library("Seurat")
#install.packages('tidyverse')
library("tidyverse")
library("Matrix")
library("scales")
library("cowplot")
library("RCurl")
options(stringsAsFactors = F)
library("scRNAseq")
library("scater")
# devtools::install_github(repo = "satijalab/seurat", ref = "loom")
library("loomR")
library("patchwork")
library("scuttle")
library("dplyr")
library("tibble")
library("HCAData")
library("SingleR")
library("org.Hs.eg.db")
library("clusterProfiler")
library("vroom")
library("celldex")
library("dittoSeq")
set.seed(2020)
library("monocle")
suppressPackageStartupMessages(library("SingleCellExperiment"))
library("ggplot2"); theme_set(theme_bw())
library("DuoClustering2018")
require(scry)
load("data/gliomaP.combined.RData")
load("data/gliomaF.combined.RData")
load("data/merged_filtered_seurat.RData")
load("data/seurat_filtered.RData")
load("data/filtered_gliomaP.RData")
load("data/filtered_gliomaF.RData")
load("data/filtered_gliomaF_norm.RData")
load("data/filtered_gliomaP_norm.RData")
load("data/filtered_seurat20210417.RData")
load("data/integrated.anchors_seurat20210417.RData")
load("data/integrated.combined_seurat20210417.RData")
load("data/glioma.combined.pca17.res0.2.RData")
load("data/glioma.pca17.markers.beforeanno.RData")
##############################################################################################################################
top30<-glioma.pca29.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_scRNA-seq/annotation/Seurat/top30.pca29.markers.csv",sep=",",quote=F)
head(Idents(glioma.combined.pca29), 5)
write.table(glioma.pca29.markers,"2_scRNA-seq/annotation/Seurat/glioma.pca29.markers.csv",sep=",",quote=F)
# Run non-linear dimensional reduction (UMAP/tSNE)
glioma.combined.pca29 <- RunUMAP(glioma.combined.pca29, dims = 1:29)
glioma.combined.pca29 <- RunTSNE(glioma.combined.pca29, dims = 1:29)
glioma.combined.pca29 <- FindNeighbors(glioma.combined.pca29, reduction = "pca", dims = 1:29)
levels(glioma.combined.pca29)
library(tidyverse)
umap_tx = glioma.combined.pca29@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(label = glioma.combined.pca29@meta.data$label) %>% cbind(subpop = glioma.combined.pca29@meta.data$integrated_snn_res.0.75)  %>% cbind(indiv = glioma.combined.pca29@meta.data$seq_folder)
write.csv(umap_tx,file="2_scRNA-seq/annotation/Seurat/umap_tx.csv",quote=F)
tsne_tx = glioma.combined.pca29@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(label = glioma.combined.pca29@meta.data$label) %>% cbind(subpop = glioma.combined.pca29@meta.data$integrated_snn_res.0.75)  %>% cbind(indiv = glioma.combined.pca29@meta.data$seq_folder)
write.csv(tsne_tx,file="2_scRNA-seq/annotation/Seurat/tsne_tx.csv",quote=F)
#umap_tx_ind = glioma.combined.pca17@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(tx = glioma.combined.pca17@meta.data$seq_folder)
#write.csv(umap_tx_ind,file="annotation/Seurat/umap_tx_ind.csv",quote=F)
#ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point() + 
#scale_color_manual(values=c("group1_untreated" = "darkblue", 
#                            "group1_treated" = "darkred"))
#tsne_tx = glioma.combined.pca17@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = glioma.combined.pca17@meta.data$label)
#tsne_tx_ind = glioma.combined.pca17@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = glioma.combined.pca17@meta.data$seq_folder)

# Visualization
pdf(paste0("2_scRNA-seq/annotation/Seurat/umap.pca29.res",r,".splitbyLabelIndiv.pdf"),width=60,height=10)
DimPlot(glioma.combined.pca29, reduction = "umap", label = TRUE, pt.size=1,label.size = 8, split.by = 'seq_folder', group.by = 'integrated_snn_res.0.75')
dev.off()
pdf(paste0("2_scRNA-seq/annotation/Seurat/umap.pca29.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(glioma.combined.pca29, reduction = "umap", label = TRUE, pt.size=1,label.size = 8, split.by = 'label', group.by = 'integrated_snn_res.0.75')
dev.off()
pdf(paste0("2_scRNA-seq/annotation/Seurat/umap.pca29.res",r,".pdf"),width=10,height=10)
DimPlot(glioma.combined.pca29, reduction = "umap", label = TRUE, pt.size=1,label.size = 8, group.by = 'integrated_snn_res.0.75')
dev.off()
pdf(paste0("2_scRNA-seq/annotation/Seurat/tsne.pca29.res",r,".splitbyLabelIndiv.pdf"),width=60,height=10)
DimPlot(glioma.combined.pca29, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, split.by = 'seq_folder', group.by = 'integrated_snn_res.0.75')
dev.off()
pdf(paste0("2_scRNA-seq/annotation/Seurat/tsne.pca29.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(glioma.combined.pca29, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, split.by = 'label', group.by = 'integrated_snn_res.0.75')
dev.off()
pdf(paste0("2_scRNA-seq/annotation/Seurat/tsne.pca29.res",r,".pdf"),width=10,height=10)
DimPlot(glioma.combined.pca29, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, group.by = 'integrated_snn_res.0.75')
dev.off()

prop.table(table(Idents(glioma.combined.pca29), glioma.combined.pca29$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca29), glioma.combined.pca29$label)))
prop.table(table(Idents(glioma.combined.pca29)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca29))))
write.csv(x = allsampleprop.each,file = '2_scRNA-seq/annotation/Seurat/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_scRNA-seq/annotation/Seurat/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(glioma.combined.pca29))
pro.total <- table(Idents(glioma.combined.pca29),glioma.combined.pca29$label)
table(Idents(glioma.combined.pca29),glioma.combined.pca29$label)
pro.each <- table(Idents(glioma.combined.pca29),glioma.combined.pca29$label)
write.csv(x =pro.total,file = '2_scRNA-seq/annotation/Seurat/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_scRNA-seq/annotation/Seurat/anno.pro.each.csv',quote = T,row.names = T)




save(glioma.combined.pca17,file="data/glioma.combined.pca17.res0.4.19clusters.aftercluster.autoSeurat.nolabel.RData")
save(glioma.combined.pca17,file="data/glioma.combined.pca17.aftercluster.autoSeurat.nolabel.RData")
prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label)))
prop.table(table(Idents(glioma.combined.pca17)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17))))
write.csv(x = allsampleprop.each,file = 'annotation/Seurat/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'annotation/Seurat/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(glioma.combined.pca17))
pro.total <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
pro.each <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
write.csv(x =pro.total,file = 'annotation/Seurat/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'annotation/Seurat/anno.pro.each.csv',quote = T,row.names = T)


##############################################################################################################################
## auto annotation by SCINA marker genes
library('SCINA')
library('preprocessCore')
DefaultAssay(glioma.combined.pca17) <- "RNA"
exp=GetAssayData(glioma.combined.pca17,slot="counts")
exp_data<-as.matrix(exp)
#Log scale and quantile normalization
exp_data=log(exp_data+1)
exp[]=normalize.quantiles(exp_data)
save(exp,file="data/exp_data.RData")
signatures=preprocess.signatures("annotation/SCINA/signatures_rmnogenes.csv")
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
    convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
#View(results$cell_labels)
unique(results$cell_labels)
#View(results$probabilities)
#pdf("annotation/annotationResultvisualization.pdf")
#plotheat.SCINA(exp, results, signatures)
#dev.off()
table(results$cell_labels)
celllabels<-results$cell_labels
cellbarcodes<-colnames(results$probabilities)
predicted<-cbind(t(t(cellbarcodes)),t(t(celllabels)))
write.table(predicted,"annotation/SCINA/predicted.txt",sep="\t",quote=F,row.names=F,col.names=F)
#predicted<-read.table("annotation/SCINA/predicted.txt",sep="\t")
glioma.combined.pca17$celltype<-predicted[,2]
glioma.combined.pca17 <- RunUMAP(glioma.combined.pca17, dims = 1:17)
glioma.combined.pca17 <- RunTSNE(glioma.combined.pca17, dims = 1:17)
pdf("annotation/SCINA/umap.SCINAanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "umap",group.by ="celltype",pt.size=2)
dev.off()
pdf("annotation/SCINA/tsne.SCINAanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "tsne",group.by ="celltype",pt.size=2)
dev.off()
save(glioma.combined.pca17,file="data/glioma.combined.pca17.afteranno.autoSCINA.RData")
rownames(predicted)<-predicted[,1]
predicted<-predicted[,-1]
predictedlabel<-factor(predicted)
Idents(glioma.combined.pca17) <- predictedlabel
save(glioma.combined.pca17,file="data/glioma.combined.pca17.afteranno.autoSCINA.labeled.RData")
prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label)))
prop.table(table(Idents(glioma.combined.pca17)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17))))
write.csv(x = allsampleprop.each,file = 'annotation/SCINA/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'annotation/SCINA/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(glioma.combined.pca17))
pro.total <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
pro.each <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
write.csv(x =pro.total,file = 'annotation/SCINA/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'annotation/SCINA/anno.pro.each.csv',quote = T,row.names = T)
glioma.combined.pca17.rmunknown<-subset(x=glioma.combined.pca17,subset = celltype != "unknown")
glioma<-glioma.combined.pca17.rmunknown
save(glioma,file="data/glioma.combined.pca17.rmunknown.RData")
#show_col(hue_pal()(11))
pdf("results/cellsubpopulations/tsne.SCINAanno.rmunknown.pdf")
DimPlot(glioma, reduction = "tsne",group.by ="celltype",pt.size=1, cols = c('Astrocyte'='#F8766D','Cancer.stem.cell'='#00A6FF','Endotheliocyte'='#DB8E00',
	'Macrophage'='#FF63B6','Microglia'='#64B200','Multilymphoid.progenitor.cell'='#00BD5C','Mural.cell'='#00C1A7','Neural.progenitor.cell'='#00BFC4',
	'Oligodendrocyte'='#B385FF','Oligodendrocyte.progenitor.cell'='#EF67EB','T.cell'='#AEA200'),label=TRUE) + NoLegend()
dev.off()
pdf("results/cellsubpopulations/tsne.SCINAanno.rmunknown.2conditions.pdf")
DimPlot(glioma, reduction = "tsne",group.by ="celltype",pt.size=1, cols = c('Astrocyte'='#F8766D','Cancer.stem.cell'='#00A6FF','Endotheliocyte'='#DB8E00',
	'Macrophage'='#FF63B6','Microglia'='#64B200','Multilymphoid.progenitor.cell'='#00BD5C','Mural.cell'='#00C1A7','Neural.progenitor.cell'='#00BFC4',
	'Oligodendrocyte'='#B385FF','Oligodendrocyte.progenitor.cell'='#EF67EB','T.cell'='#AEA200'),label=TRUE,split.by = "label") + NoLegend()
dev.off()


prop.table(table(Idents(glioma), glioma$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma), glioma$label)))
prop.table(table(Idents(glioma)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma))))
write.csv(x = allsampleprop.each,file = 'results/cellsubpopulations/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'results/cellsubpopulations/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(glioma))
pro.total <- table(Idents(glioma),glioma$label)
table(Idents(glioma),glioma$label)
pro.each <- table(Idents(glioma),glioma$label)
write.csv(x =pro.total,file = 'results/cellsubpopulations/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'results/cellsubpopulations/anno.pro.each.csv',quote = T,row.names = T)
pro.total<-read.csv("results/cellsubpopulations/anno.pro.total.csv",row.names = 1)
pdf("results/cellsubpopulations/statistics.gliomaP.pdf")
ggplot(pro.total,aes(x=row.names(pro.total),y=gliomaP)) + geom_bar(stat='identity') + coord_flip() + ylim(0,25000)
dev.off()
pdf("results/cellsubpopulations/statistics.gliomaF.pdf")
ggplot(pro.total,aes(x=row.names(pro.total),y=gliomaF)) + geom_bar(stat='identity') + coord_flip() + ylim(0,25000) 
dev.off()
pro.total.sample <- table(Idents(glioma),glioma$seq_folder)
write.csv(x =pro.total.sample,file = 'results/cellsubpopulations/anno.pro.total.sample.csv',quote = T,row.names = T)
library(reshape2)
library(dplyr)
pro.total.sample<-read.csv("results/cellsubpopulations/anno.pro.total.sample.csv",row.names = 1)
pro.total.sample_rownames <- rownames(pro.total.sample)
pro.total.sample_colnames <- colnames(pro.total.sample)
pro.total.sample$type <- pro.total.sample_rownames
pro.total.sample_m <- melt(pro.total.sample, id.vars=c("type"))
pdf("results/cellsubpopulations/statistics.pdf")
p <- ggplot(pro.total.sample_m, aes(x=variable, y=value)) +
  geom_bar(stat="identity", position="fill", aes(fill=type))
p
dev.off()
library(psych)
library(pheatmap)
glioma.markers <- FindAllMarkers(object = glioma, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(glioma.markers,"results/markergenes/glioma.markers.csv",sep=",",quote=F)
top30<-glioma.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"results/markergenes/top30.markers.csv",sep=",",quote=F)
AverageExp<-AverageExpression(glioma,features=unique(top30$gene))
typeof(AverageExp)
head(AverageExp$RNA)
DefaultAssay(glioma) <- "integrated"
#levels(glioma)
#levels(glioma) <- c("Macrophage", "Astrocyte", "Oligodendrocyte", "T.cell", "Oligodendrocyte.progenitor.cell","Neural.progenitor.cell",
	"Cancer.stem.cell", "Mural.cell","Multilymphoid.progenitor.cell","Microglia","Endotheliocyte")
#levels(glioma)
pdf("results/markergenes/averageExptop30.clusters.pdf")
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r,cluster_row = FALSE,cluster_col = FALSE)
dev.off()
pdf("results/markergenes/heatmap.top30.pdf",width=24,height=18)
DoHeatmap(glioma,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()
levels(glioma)
levels(glioma) <- c("Astrocyte","Cancer.stem.cell","Endotheliocyte","Macrophage","Microglia","Multilymphoid.progenitor.cell",
	"Mural.cell","Neural.progenitor.cell","Oligodendrocyte","Oligodendrocyte.progenitor.cell","T.cell")
levels(glioma)
DefaultAssay(glioma) <- "RNA"
features.plot <- c("CLU","MT3","CCND2","SOX4","CLDN5","ITM2A","HLA-DRA","HLA-DRB1",
	"LYZ","S100A8","ZFAND2A","VEGFA","DCN","COL3A1","CENPF","TOP2A","PLP1","TF","BCAN","OLIG1","GNLY","CCL5")
pdf("results/markergenes/markergenes.dotplot.pdf",width = 10, height = 8)
DotPlot(object = glioma, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()


pdf("results/markergenes/markergenes.Macrophage.pdf")
FeaturePlot(glioma, features = "CD14", reduction="tsne")
dev.off()
pdf("results/markergenes/markergenes.Cancerstemcell.pdf")
FeaturePlot(glioma, features = "SOX2", reduction="tsne")
dev.off()

features.plot <- c("HLA-DQB1","FOLR2","PTN","EGFR","PLP1","MAG","HLA-DRB5","IL1B","RAMP1","IGFBP7","NKG7","CCL5","HSPA6","ZFAND2A","CENPF","TOP2A","DCN","COL1A2","BCAS1","FYN","CLDN5","ANGPT2")
pdf("results/pca16res0.3/test.dotplot.pdf",width = 10, height = 8)
DotPlot(object = glioma.combined.pca16, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()

##############################################################################################################################
# auto annotation by SingleR training dataset
glioma.combined.pca17@meta.data$cell.type <- Idents(glioma.combined.pca17)
test <- as.SingleCellExperiment(glioma.combined.pca17)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData()
sceESC<-LaMannoBrainData('human-es')
sceEmidBrain<-LaMannoBrainData('human-embryo')
sceIPSC<-LaMannoBrainData('human-ips')
sceESC<-sceESC[,!is.na(sceESC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceESC <- logNormCounts(sceESC)
sceEmidBrain<-sceEmidBrain[,!is.na(sceEmidBrain$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceEmidBrain <- logNormCounts(sceEmidBrain)
sceIPSC<-sceIPSC[,!is.na(sceIPSC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceIPSC <- logNormCounts(sceIPSC)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se, ESC=sceESC, midbrain=sceEmidBrain, iPSC=sceIPSC),
                labels = list(hpca.se$label.main , bpe.se$label.main, sceESC$Cell_type, sceEmidBrain$Cell_type, sceIPSC$Cell_type),
                method = "cluster",
                cluster = test$cell.type)
table(Anno$pruned.labels)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
#将细胞注释信息重新添加到Seurat对象中去
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(glioma.combined.pca17)
glioma.combined.pca17 <- RenameIdents(glioma.combined.pca17, new.cluster.ids)
head(Idents(glioma.combined.pca17), 5)
levels(Idents(glioma.combined.pca17))
glioma.combined.pca17$celltype<-Idents(glioma.combined.pca17)
# Run non-linear dimensional reduction (UMAP/tSNE)
glioma.combined.pca17 <- RunUMAP(glioma.combined.pca17, reduction = "pca", dims = 1:16)
save(glioma.combined.pca17,file="data/glioma.combined.pca17.res0.3.afteranno.autoSingleR.RData")
# Visualization
glioma.combined.pca17 <- RunUMAP(glioma.combined.pca17, dims = 1:17)
glioma.combined.pca17 <- RunTSNE(glioma.combined.pca17, dims = 1:17)
pdf("annotation/SingleR/umap.singleRanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "umap",group.by ="celltype",pt.size=2)
dev.off()
pdf("annotation/SingleR/tsne.singleRanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "tsne",group.by ="celltype",pt.size=2)
dev.off()
prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label)))
prop.table(table(Idents(glioma.combined.pca17)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17))))
write.csv(x = allsampleprop.each,file = 'annotation/SingleR/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'annotation/SingleR/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(glioma.combined.pca17))
pro.total <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
pro.each <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
write.csv(x =pro.total,file = 'annotation/SingleR/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'annotation/SingleR/anno.pro.each.csv',quote = T,row.names = T)
##############################################################################################################################
glioma.pca16.markers <- FindAllMarkers(object = glioma.combined.pca16, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top30<-glioma.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_logFC)
write.table(top30,"results/pca16res0.3/glioma.pca16.top30.anno.csv",sep=",",quote=F)
glioma.combined.pca16$celltype<-Idents(glioma.combined.pca16)
pdf("results/pca16res0.3/heatmap.top30.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca16,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()



new.cluster.ids<-c("Macrophage", "Cancer stem cell", "Macrophage", "Oligodendrocyte", "Microglia", "Astrocyte", "Macrophage", "T cell", "Multilymphoid progenitor cell", "Neural progenitor cell", "Mural cell", "Neural progenitor cell", "Oligodendrocyte progenitor cell", "Endotheliocyte")
names(new.cluster.ids) <- levels(glioma.combined.pca16)
glioma.combined.pca16 <- RenameIdents(glioma.combined.pca16, new.cluster.ids)
glioma.combined.pca16 <- RunUMAP(glioma.combined.pca16, dims = 1:16)
save(glioma.combined.pca16,file="data/glioma.combined.pca16.res0.3.afteranno.RData")
glioma.pca16.markers <- FindAllMarkers(object = glioma.combined.pca16, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(glioma.pca16.markers,"results/pca16res0.3/glioma.pca16.markers.anno.csv",sep=",",quote=F)
save(glioma.pca16.markers,file="data/glioma.pca16.markers.afteranno.RData")

pdf("results/pca16res0.3/umap.pca16.res0.3.integrate.pdf", width = 36, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label",pt.size=2)
p2<-DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE,pt.size=2)
plot_grid(p1, p2)
dev.off()
pdf("results/pca16res0.3/umap.pca16.cluster.res0.3.2treat.pdf", width = 36, height = 18)
DimPlot(glioma.combined.pca16, reduction = "umap", split.by = "label",pt.size=2)
dev.off()
top30<-glioma.pca16.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_logFC)
write.table(top30,"results/pca16res0.3/glioma.pca16.top30.anno.csv",sep=",",quote=F)
glioma.combined.pca16$celltype<-Idents(glioma.combined.pca16)
pdf("results/pca16res0.3/heatmap.top30.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca16,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()

prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label)))
prop.table(table(Idents(glioma.combined.pca16)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16))))
write.csv(x = allsampleprop.each,file = 'results/pca16res0.3/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'results/pca16res0.3/anno.allsample.total.prop.csv',quote = T,row.names = T)


table(Idents(glioma.combined.pca16))
pro.total <- table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
pro.each <- table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
write.csv(x =pro.total,file = 'results/pca16res0.3/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'results/pca16res0.3/anno.pro.each.csv',quote = T,row.names = T)


delta.genes <- c("GPNMB","PLIN2","EGFR","MEG3","PLP1","PPP1R14A","HLA-DRB5","IL1B","CHI3L1","CXCL14","NKG7","GNLY","IGF1","EGR3","HIST1H4C","CENPF","DCN","COL3A1","BCAS1","RP11-89N17.4","CLDN5","ITM2A")
pdf("results/pca16res0.3/dittoDotPlot.pdf")
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()


delta.genes <- c("SPP1","GPNMB","APOC1","FTL","CTSB","S100A9","HMOX1","CD14","HLA-DRA","CTSD","C1QB","C1QA","HLA-DPA1","CTSL","ALOX5AP","LGMN","S100A8","HLA-DQB1","NUPR1","FOLR2","HLA-DPB1","HLA-DRB1","MS4A4A","MT1G","HLA-DQA1","STAB1","SERPINA1","TGFBI","PLIN2","RNASE1")
pdf("results/pca16res0.3/dittoDotPlot.Macrophage.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("SEC61G","PTN","EGFR","PTPRZ1","C1orf61","SOX4","MEG3","NOVA1","TUBB2B","BCAN","SOX2","IGFBP2","C1QL1","GRIA2","CKB","CLU","FABP7","GPM6A","CCND2","FAM181B","SERPINE2","OLIG1","METRN","PCSK1N","BEX1","LINC00461","SCRG1","VGF","ID4","VEGFA")
pdf("results/pca16res0.3/dittoDotPlot.Cancerstemcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("PLP1","TF","PPP1R14A","PTGDS","MBP","NKX6-2","MAG","QDPR","CNP","CLDND1","CRYAB","APLP1","ERMN","CLDN11","CNDP1","SEPP1","KLK6","TMEM144","GPM6B","MOG","ABCA2","SLC44A1","CAMK2N1","AMER2","EDIL3","S100A1","ENPP2","KCNMB4","ANLN","SLAIN1")
pdf("results/pca16res0.3/dittoDotPlot.Oligodendrocyte.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("HLA-DRB5","CCL3","CCL4","IL1B","CCL3L1","HLA-DRB1","IER3","CD83","HLA-DPB1","IL8","FCGBP","NR4A2","CCL3L3","HLA-DPA1","HLA-DRA","C1QA","HLA-DQA1","CD163","C1QB","CXCL3","HLA-DQB1","TNF","EGR2","NR4A1","CCL4L2","CCL2","ZNF331","LYZ","CXCL2","CD14")
pdf("results/pca16res0.3/dittoDotPlot.Microglia.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("CHI3L1","MT3","MT2A","GFAP","RAMP1","IGFBP7","HILPDA","CLU","IGFBP5","S100A10","RARRES2","MGST1","CXCL14","CSRP2","VIM","CPE","S100A16","MAP1B","SPARCL1","ANXA1","IGFBP2","CAV1","MT1X","DBI","AKAP12","MRPS6","OCIAD2","LGALS3","GAP43","IGFBP3")
pdf("results/pca16res0.3/dittoDotPlot.Astrocyte.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("NKG7","CCL5","IL32","GZMA","CD2","CD52","GZMH","CD3D","GZMB","PTPRCAP","CST7","CD3G","PRF1","RAC2","CD3E","ETS1","S100A4","CD7","TXNIP","DUSP2","RUNX3","IL2RG","CD69","GZMK","IL7R","SYNE2","CD8A","GNLY","KLRD1","FGFBP2")
pdf("results/pca16res0.3/dittoDotPlot.Tcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("HSPA6","ZFAND2A","FAM46A","VEGFA","DDIT3","IFI30","SPRY1","TNF","OTUD1","MDM2","CDKN1A","IGF1","PDK4","HSPH1","KMT2E-AS1","EGR3","CCL4L2","HSPA1A","SNHG12","KLF2","GADD45G","PMAIP1","HSPB1","HSPA1B","RGS16","DNAJB1","NRP2","RGS2","PLK2","CYB5D2")
pdf("results/pca16res0.3/dittoDotPlot.Multilymphoidprogenitorcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("HIST1H4C","CENPF","TOP2A","MKI67","TUBB","HMGB2","HMGN2","H2AFZ","UBE2S","UBE2C","TUBA1B","BIRC5","NUSAP1","ASPM","PTTG1","TYMS","KIAA0101","TPX2","STMN1","HIST1H1B","PRC1","MAD2L1","CDK1","CENPE","SEC61G","KPNA2","CCNB1","CDKN3","SMC4","VGF")
pdf("results/pca16res0.3/dittoDotPlot.Neuralprogenitorcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("DCN","COL3A1","COL1A1","APOD","COL1A2","COL6A3","CALD1","LUM","COL6A2","MGP","IGFBP7","NDUFA4L2","COL4A1","IGF2","COL6A1","FN1","CAV1","TIMP3","COL4A2","PCOLCE","NR2F2","LAMB1","CFH","NID1","PDGFRB","PTGDS","C1R","RARRES2","S100A6","TIMP1")
pdf("results/pca16res0.3/dittoDotPlot.Muralcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("BCAS1","FXYD6","SCRG1","FYN","CDKN1C","SIRT2","BCAN","RP11-161M6.2","MARCKSL1","SMOC1","SOX8","TNR","GPM6A","SOX4","DLL3","SHD","GRIA2","NKAIN4","MDFI","GPR17","RND2","FERMT1","NCAM1","SEC61G","BEX1","SOX6","ATCAY","OLIG1","LSAMP","RP11-89N17.4")
pdf("results/pca16res0.3/dittoDotPlot.Oligodendrocyteprogenitorcell.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

delta.genes <- c("PODXL","CLDN5","GPR116","SLC9A3R2","SPARCL1","SLC7A5","EPAS1","ESAM","TIMP3","SDPR","ITM2A","IGFBP7","EGFL7","INSR","FLT1","IFI27","TSC22D1","SLC38A5","CTGF","VWF","ANGPT2","MFSD2A","ABCG2","SPARC","GNG11","COL4A1","BSG","RGCC","EDN1","ADIRF")
pdf("results/pca16res0.3/dittoDotPlot.Endotheliocyte.pdf", width = 54, height = 18)
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()



delta.genes <- c("HLA-DQB1","FOLR2","PTN","EGFR","PLP1","MAG","HLA-DRB5","IL1B","RAMP1","IGFBP7","NKG7","CCL5","HSPA6","ZFAND2A","CENPF","TOP2A","DCN","COL1A2","BCAS1","FYN","CLDN5","ANGPT2")
pdf("results/pca16res0.3/dittoDotPlot.revised.pdf")
dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()

features.plot <- c("HLA-DQB1","FOLR2","PTN","EGFR","PLP1","MAG","HLA-DRB5","IL1B","RAMP1","IGFBP7","NKG7","CCL5","HSPA6","ZFAND2A","CENPF","TOP2A","DCN","COL1A2","BCAS1","FYN","CLDN5","ANGPT2")
pdf("results/pca16res0.3/test.dotplot.pdf",width = 10, height = 8)
DotPlot(object = glioma.combined.pca16, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()

features.plot <- c("HLA-DQB1","PTN","PLP1","HLA-DRB5","RAMP1","NKG7","HSPA6","CENPF","DCN","BCAS1","CLDN5")
pdf("results/pca16res0.3/test.RidgePlot.pdf",width = 20, height = 20)
RidgePlot(object = glioma.combined.pca16, features = features.plot,ncol=3)
dev.off()



features.plot <- c("HLA-DQB1","FOLR2","MS4A4A","HLA-DQA1","STAB1","SERPINA1")
pdf("results/pca16res0.3/macrophage.RidgePlot.pdf",width = 20, height = 10)
RidgePlot(object = glioma.combined.pca16, features = features.plot,ncol=3)
dev.off()
features.plot <- c("PTN","EGFR","MEG3","SOX2","CCND2","BEX1")
pdf("results/pca16res0.3/cancerstemcell.RidgePlot.pdf",width = 20, height = 10)
RidgePlot(object = glioma.combined.pca16, features = features.plot,ncol=3)
dev.off()




dittoDotPlot(glioma.combined.pca16, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()



prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label)))
prop.table(table(Idents(glioma.combined.pca16)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16))))
write.csv(x = allsampleprop.each,file = 'results/pca16res0.3/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'results/pca16res0.3/anno.allsample.total.prop.csv',quote = T,row.names = T)

table(Idents(glioma.combined.pca16))
pro.total <- table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
pro.each <- table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
write.csv(x =pro.total,file = 'results/pca16res0.3/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'results/pca16res0.3/anno.pro.each.csv',quote = T,row.names = T)

AverageExp<-AverageExpression(glioma.combined.pca16,features=unique(top30$gene))
typeof(AverageExp)
head(AverageExp$RNA)
library(psych)
library(pheatmap)
pdf("results/pca16res0.3/averageExptop30.clusters.pdf")
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)
dev.off()


require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(glioma.combined.pca16))
names(geneList)<-row.names(glioma.combined.pca16)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(glioma.pca16.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(glioma.pca16.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="data/gene_erichment_results.RData")
#Macrophage
gene_erichment_results[["Macrophage"]][["topGO"]][1:5,]
pdf("results/pca16res0.3/clusterAnalysis/Macrophage.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["Macrophage"]][["DGN"]], showCategory=30) 
dev.off()
## categorySize can be scaled by 'pvalue' or 'geneNum'
pdf("results/pca16res0.3/clusterAnalysis/Macrophage.2.pdf")
p1<-cnetplot(gene_erichment_results[["Macrophage"]][["DGN"]], categorySize="pvalue", foldChange=geneList)
p2<-cnetplot(gene_erichment_results[["Macrophage"]][["DGN"]], foldChange=geneList, circular = TRUE, colorEdge = TRUE)
plot_grid(p1, p2, ncol=2)
dev.off()
#Tcell
gene_erichment_results[["T cell"]][["topGO"]][1:5,]
pdf("results/pca16res0.3/clusterAnalysis/T cell.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["T cell"]][["DGN"]], showCategory=30) 
dev.off()
#Cancer stem cell
gene_erichment_results[["Cancer stem cell"]][["topGO"]][1:5,]
pdf("results/pca16res0.3/clusterAnalysis/Cancer stem cell.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["Cancer stem cell"]][["DGN"]], showCategory=30) 
dev.off()
#Microglia
gene_erichment_results[["Microglia"]][["topGO"]][1:5,]
pdf("results/pca16res0.3/clusterAnalysis/Microglia.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["Microglia"]][["DGN"]], showCategory=30) 
dev.off()


#Macrophage
pdf("results/pca16res0.3/clusterAnalysis/markergenes.Macrophage.pdf",width=15,height=10)
VlnPlot(glioma.combined.pca16, features = c("GPNMB","PLIN2"), group.by = 'celltype', pt.size = 0) 
dev.off()
#https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
FeaturePlot(glioma.combined.pca16, features = c("CD14","Sox2"))

#Tcell
pdf("results/pca16res0.3/clusterAnalysis/markergenes.Macrophage.pdf",width=15,height=10)
VlnPlot(glioma.combined.pca16, features = c("GPNMB","PLIN2"), group.by = 'celltype', pt.size = 0) 
dev.off()
#https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
FeaturePlot(glioma.combined.pca16, features = c("CD14","Sox2"))

#Cancer stem cell
pdf("results/pca16res0.3/clusterAnalysis/markergenes.Cancerstemcell.pdf")
VlnPlot(glioma.combined.pca16, features = c("EGFR","MEG3"), group.by = 'celltype', pt.size = 0) 
dev.off()
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737740/

pdf("results/pca16res0.3/markergenes/markergenes.distribution.pdf",width=20,height=10)
FeaturePlot(glioma.combined.pca16, features = c("CD14","SOX2"))
dev.off()


pdf("TF/SPI1.RELA.distribution.pdf",width=20,height=10)
FeaturePlot(glioma.combined.pca16, features = c("SPI1","RELA"))
dev.off()

pdf("TF/CEBPA.CEBPB.distribution.pdf",width=20,height=10)
FeaturePlot(glioma.combined.pca16, features = c("CEBPA","CEBPB"))
dev.off()

pdf("TF/MAX.NR4A1.distribution.pdf",width=20,height=10)
FeaturePlot(glioma.combined.pca16, features = c("MAX","NR4A1"))
dev.off()


pdf("results/pca16res0.3/test.pdf")
DoHeatmap(glioma.combined.pca16, features = VariableFeatures(glioma.combined.pca16)[1:10], cells = 1:500, size = 4, angle = 90) +  NoLegend()
dev.off()



pdf("figures/heatmap.top10.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca18,features=top10$gene,cells = 1:500, size = 4, angle = 90) + NoLegend()
dev.off()
save(glioma.pca18.markers,file="data/glioma.pca18.markers.namedcelltypes.RData")
save(glioma.combined.pca18,file="data/glioma.combined.pca18.namedcelltypes.RData")

top10<-glioma.pca18.markers %>% group_by(seurat_clusters) %>% top_n(n=10,wt=avg_logFC)
pdf("figures/heatmap.top10.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca18,features=top10$gene,cells = 1:500, size = 4, angle = 90) + NoLegend()
dev.off()


#https://www.bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html

if(FALSE){
glioma.combined.pca17@meta.data$cell.type <- Idents(glioma.combined.pca17)
test <- as.SingleCellExperiment(glioma.combined.pca17)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData()
sceESC<-LaMannoBrainData('human-es')
sceEmidBrain<-LaMannoBrainData('human-embryo')
sceIPSC<-LaMannoBrainData('human-ips')
sceESC<-sceESC[,!is.na(sceESC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceESC <- logNormCounts(sceESC)
sceEmidBrain<-sceEmidBrain[,!is.na(sceEmidBrain$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceEmidBrain <- logNormCounts(sceEmidBrain)
sceIPSC<-sceIPSC[,!is.na(sceIPSC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceIPSC <- logNormCounts(sceIPSC)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se, ESC=sceESC, midbrain=sceEmidBrain, iPSC=sceIPSC),
                labels = list(hpca.se$label.main , bpe.se$label.main, sceESC$Cell_type, sceEmidBrain$Cell_type, sceIPSC$Cell_type),
                method = "cluster",
                cluster = test$cell.type)
table(Anno$pruned.labels)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
#将细胞注释信息重新添加到Seurat对象中去
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(glioma.combined.pca17)
glioma.combined.pca17 <- RenameIdents(glioma.combined.pca17, new.cluster.ids)
head(Idents(glioma.combined.pca17), 5)
levels(Idents(glioma.combined.pca17))
glioma.combined.pca17$celltype<-Idents(glioma.combined.pca17)
# Run non-linear dimensional reduction (UMAP/tSNE)
glioma.combined.pca17 <- RunUMAP(glioma.combined.pca17, reduction = "pca", dims = 1:16)
save(glioma.combined.pca17,file="data/glioma.combined.pca17.res0.3.afteranno.autoSingleR.RData")
# Visualization

glioma.combined.pca17 <- RunUMAP(glioma.combined.pca17, dims = 1:17)
glioma.combined.pca17 <- RunTSNE(glioma.combined.pca17, dims = 1:17)

pdf("umap.singleRanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "umap",group.by ="celltype",pt.size=2)
dev.off()
pdf("tsne.singleRanno.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "tsne",group.by ="celltype",pt.size=2)
dev.off()



pdf("umap.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "umap",pt.size=2)
dev.off()
pdf("tsne.pdf",width=10,height=10)
DimPlot(glioma.combined.pca17, reduction = "tsne",pt.size=2)
dev.off()

prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17), glioma.combined.pca17$label)))
prop.table(table(Idents(glioma.combined.pca17)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca17))))
write.csv(x = allsampleprop.each,file = 'results/pca16res0.3/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'results/pca16res0.3/anno.allsample.total.prop.csv',quote = T,row.names = T)



table(Idents(glioma.combined.pca17))
pro.total <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
pro.each <- table(Idents(glioma.combined.pca17),glioma.combined.pca17$label)
write.csv(x =pro.total,file = 'results/pca16res0.3/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = 'results/pca16res0.3/anno.pro.each.csv',quote = T,row.names = T)

pdf("results/pca16res0.3/UMAP.singleRanno.pdf", width = 24, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label",pt.size = 2)
p2 <- DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE,pt.size = 2)
plot_grid(p1, p2)
dev.off()
pdf(paste0("figures/umap.pca16.res",r,".cluster.pdf"), width = 24, height = 18)
DimPlot(glioma.combined.pca16, reduction = "umap", split.by = "label")
dev.off()

glioma.pca16.markers <- FindAllMarkers(object = glioma.combined.pca16, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
save(glioma.pca18.markers,file="data/glioma.pca18.markers.beforecelltypes.RData")

top10<-glioma.pca18.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
pdf("figures/heatmap.top10.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca18,features=top10$gene,cells = 1:500, size = 4, angle = 90) + NoLegend()
dev.off()
save(glioma.pca18.markers,file="data/glioma.pca18.markers.namedcelltypes.RData")
save(glioma.combined.pca18,file="data/glioma.combined.pca18.namedcelltypes.RData")

top10<-glioma.pca18.markers %>% group_by(seurat_clusters) %>% top_n(n=10,wt=avg_logFC)
pdf("figures/heatmap.top10.pdf",width=24,height=18)
DoHeatmap(glioma.combined.pca18,features=top10$gene,cells = 1:500, size = 4, angle = 90) + NoLegend()
dev.off()

head(glioma.pca18.markers)
gene.df <- bitr(rownames(glioma.pca18.markers), fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的  ENTREZID
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb =org.Hs.eg.db)
head(gene.df)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers
y <- enricher(gene.df$ENTREZID, TERM2GENE=cell_markers, minGSSize=1)
#DT::datatable(as.data.frame(y))
id <- gene.df$ENTREZID
names(id) <- gene.df$SYMBOL  # 用R语言的names函数构建类似字典的数据结构
id
y@result-> res  # 以防数据被搞坏，新建一个
res$genesym <- unlist(lapply(y@result$geneID,  FUN =function(x){paste( unlist(lapply(unlist(str_split(x,"/")),FUN=function(x){names(id[which(id ==x)])})) , collapse = "/")} ))
#DT::datatable(res)
id <- unique(gene.df$ENTREZID)
names(id) <- unique(gene.df$SYMBOL)
id
y@result-> res
res$genesym <- unlist(lapply(res$geneID,  FUN =function(x){paste( unlist(lapply(unlist(str_split(x,"/")),FUN=function(x){names(id[which(id ==x)])})) , collapse = "/")} ))
#DT::datatable(res)
write.table(res,"data/annotation.csv",sep=",")
}

glioma.combined.pca16 <- RunUMAP(glioma.combined.pca16, reduction = "pca", dims = 1:16)
glioma.combined.pca16 <- FindNeighbors(glioma.combined, reduction = "pca", dims = 1:16)
r=0.3
glioma.combined.pca16 <- FindClusters(glioma.combined.pca16, resolution = r)
levels(glioma.combined.pca16)
pdf("results/pca16res0.3/figure1.pca16.res0.3.cluster14.pdf", width = 30, height = 20)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label")
p2 <- DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
pdf("results/pca16res0.3/figure2.pca16.res0.3.cluster14.pdf", width = 30, height = 20)
DimPlot(glioma.combined.pca16, reduction = "umap", split.by = "label")
dev.off()
glioma.pca18.markers <- FindAllMarkers(object = glioma.combined.pca18, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
head(Idents(glioma.combined.pca18), 5)
glioma.combined.pca18@meta.data$cell.type <- Idents(glioma.combined.pca18)
save(glioma.combined.pca18,file="data/glioma.combined.pca18.RData")
save(glioma.pca18.markers,file="data/glioma.pca18.markers.RData")

data <- as(as.matrix(glioma.combined.pca16@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = glioma.combined.pca16@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
monocle_cds@phenoData@data[["cell_type"]] <- as.character(glioma.combined.pca16@active.ident)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
#Filtering low-quality cells
monocle_cds <- detectGenes(monocle_cds, min_expr = 3 )
print(head(fData(monocle_cds)))
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))
print(head(pData(monocle_cds)))
save(monocle_cds,file="data/monocle_cds.RData")
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],fullModelFormulaStr = "~cell_type")
head(diff_test_res)
write.table(diff_test_res,"results/pca16res0.3/pseudotime/diff_test_res.txt",sep="\t",quote=F)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 
length(ordering_genes)
write.table(ordering_genes,"results/pca16res0.3/pseudotime/ordering_genes.txt",sep="\t",quote=F,row.names=F,col.names=F)
dev <- setOrderingFilter(monocle_cds, ordering_genes)
pdf("results/pca16res0.3/pseudotime/ordering_gene.pdf", width = 24, height = 18)
plot_ordering_genes(dev)
dev.off()
#Trajectory step 2: reduce data dimensionality  
dev <- reduceDimension(dev, max_components = 2,method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory  
dev <- orderCells(dev)
pdf("results/pca16res0.3/pseudotime/plot_cell_trajectory.pdf", width = 24, height = 18)
plot_cell_trajectory(gliomaP_dev, color_by = "seurat_clusters")
plot_cell_trajectory(gliomaP_dev, color_by = "Pseudotime")
plot_cell_trajectory(gliomaP_dev,color_by="cell_type")
plot_cell_trajectory(gliomaP_dev,color_by="label")+facet_wrap(~label, nrow = 1)
dev.off()

glioma.combined.pca16_gliomaP<-subset(x=glioma.combined.pca16,subset = label == "gliomaP")
glioma.combined.pca16_gliomaF<-subset(x=glioma.combined.pca16,subset = label == "gliomaF")

data <- as(as.matrix(glioma.combined.pca16_gliomaP@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = glioma.combined.pca16_gliomaP@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_gliomaP_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
monocle_gliomaP_cds@phenoData@data[["cell_type"]] <- as.character(glioma.combined.pca16_gliomaP@active.ident)
monocle_gliomaP_cds <- estimateSizeFactors(monocle_gliomaP_cds)
monocle_gliomaP_cds <- estimateDispersions(monocle_gliomaP_cds)
#Filtering low-quality cells
monocle_gliomaP_cds <- detectGenes(monocle_gliomaP_cds, min_expr = 3 )
print(head(fData(monocle_gliomaP_cds)))
expressed_genes <- row.names(subset(fData(monocle_gliomaP_cds),num_cells_expressed >= 10))
print(head(pData(monocle_gliomaP_cds)))
save(monocle_gliomaP_cds,file="data/monocle_gliomaP_cds.RData")
diff_test_res_gliomaP <- differentialGeneTest(monocle_gliomaP_cds[expressed_genes,],fullModelFormulaStr = "~cell_type")
head(diff_test_res_gliomaP)
write.table(diff_test_res_gliomaP,"results/pca16res0.3/pseudotime/diff_test_res_gliomaP.txt",sep="\t",quote=F)
ordering_genes_gliomaP <- row.names (subset(diff_test_res_gliomaP, qval < 0.01)) 
length(ordering_genes_gliomaP)
write.table(ordering_genes_gliomaP,"results/pca16res0.3/pseudotime/ordering_genes_gliomaP.txt",sep="\t",quote=F,row.names=F,col.names=F)
gliomaP_dev <- setOrderingFilter(monocle_gliomaP_cds, ordering_genes_gliomaP)
pdf("results/pca16res0.3/pseudotime/ordering_gene_gliomaP.pdf", width = 24, height = 18)
plot_ordering_genes(gliomaP_dev)
dev.off()
#Trajectory step 2: reduce data dimensionality  
gliomaP_dev <- reduceDimension(gliomaP_dev, max_components = 2,method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory  
gliomaP_dev <- orderCells(gliomaP_dev)
save(gliomaP_dev,file="data/gliomaP_dev.RData")
pdf("results/pca16res0.3/pseudotime/plot_cell_trajectory_gliomaP.pdf", width = 24, height = 18)
plot_cell_trajectory(gliomaP_dev, color_by = "seurat_clusters")
plot_cell_trajectory(gliomaP_dev, color_by = "Pseudotime")
plot_cell_trajectory(gliomaP_dev,color_by="cell_type")
plot_cell_trajectory(gliomaP_dev,color_by="label")+facet_wrap(~label, nrow = 1)
dev.off()
pdf("results/pca16res0.3/pseudotime/genes.discelltype.gliomaP.pdf", width = 27, height = 18)
to_be_tested <- row.names(subset(fData(monocle_gliomaP_cds),
              gene_short_name %in% c("APOC1", "CD14", "GPNMB", "PLIN2")))
cds_subset_gliomaP <- monocle_gliomaP_cds[to_be_tested,]
plot_genes_jitter(cds_subset_gliomaP,
                  grouping = "cell_type",
                  color_by = "cell_type",
                  nrow= 2,
                  ncol = NULL,
                  plot_trend = TRUE)
dev.off()
pdf("results/pca16res0.3/pseudotime/genes.changeasafunctionofpsedotime.gliomaP.pdf", width = 27, height = 18)
diff_test_res <- differentialGeneTest(gliomaP_dev,fullModelFormulaStr = "~sm.ns(Pseudotime)")
save(diff_test_res,file="data/diff_test_res_gliomaP.Rdata")
diff_test_res[,c("gene_short_name", "pval", "qval")]
diff_test_res_filterbyExp <- subset(diff_test_res, num_cells_expressed >= 10)
#diff_test_res_filterbySig <- subset(diff_test_res_filterbyExp, qval < 10e-10)
sig_gene_names <- row.names(subset(diff_test_res_filterbyExp, qval < 10e-10))
pdf("results/pca16res0.3/pseudotime/heatmap.gliomaP.pdf")
plot_pseudotime_heatmap(gliomaP_dev[sig_gene_names,],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)
dev.off()
                      


BEAM_res <- BEAM(gliomaP_dev, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf("results/pca16res0.3/pseudotime/plot_genes_branched_heatmap.gliomaP.pdf")
plot_genes_branched_heatmap(gliomaP_dev[row.names(subset(BEAM_res,
                                          qval < 1e-156)),],
                                          branch_point = 1,
                                          num_clusters = 3,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T)                      
dev.off()                      
                        
to_be_tested <- row.names(subset(fData(monocle_gliomaP_cds),
              gene_short_name %in% c("APOC1", "CD14", "GPNMB", "PLIN2")))
cds_subset_gliomaP <- monocle_gliomaP_cds[to_be_tested,]
plot_genes_jitter(cds_subset_gliomaP,
                  grouping = "cell_type",
                  color_by = "cell_type",
                  nrow= 2,
                  ncol = NULL,
                  plot_trend = TRUE)
dev.off()


data <- as(as.matrix(glioma.combined.pca16_gliomaF@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = glioma.combined.pca16_gliomaF@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_gliomaF_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
monocle_gliomaF_cds@phenoData@data[["cell_type"]] <- as.character(glioma.combined.pca16_gliomaF@active.ident)
monocle_gliomaF_cds <- estimateSizeFactors(monocle_gliomaF_cds)
monocle_gliomaF_cds <- estimateDispersions(monocle_gliomaF_cds)
#Filtering low-quality cells
monocle_gliomaF_cds <- detectGenes(monocle_gliomaF_cds, min_expr = 3 )
print(head(fData(monocle_gliomaF_cds)))
expressed_genes <- row.names(subset(fData(monocle_gliomaF_cds),num_cells_expressed >= 10))
print(head(pData(monocle_gliomaF_cds)))
save(monocle_gliomaF_cds,file="data/monocle_gliomaF_cds.RData")
diff_test_res_gliomaF <- differentialGeneTest(monocle_gliomaF_cds[expressed_genes,],fullModelFormulaStr = "~cell_type")
head(diff_test_res_gliomaF)
write.table(diff_test_res_gliomaF,"results/pca16res0.3/pseudotime/diff_test_res_gliomaF.txt",sep="\t",quote=F)
ordering_genes_gliomaF <- row.names (subset(diff_test_res_gliomaF, qval < 0.01)) 
ordering_genes_gliomaF <- row.names (subset(diff_test_res_gliomaF, qval < 0.01)) 
length(ordering_genes_gliomaF)
write.table(ordering_genes_gliomaF,"results/pca16res0.3/pseudotime/ordering_genes_gliomaF.txt",sep="\t",quote=F,row.names=F,col.names=F)
gliomaF_dev <- setOrderingFilter(monocle_gliomaF_cds, ordering_genes_gliomaF)
pdf("results/pca16res0.3/pseudotime/ordering_gene_gliomaF.pdf", width = 24, height = 18)
plot_ordering_genes(gliomaF_dev)
dev.off()
#Trajectory step 2: reduce data dimensionality  
gliomaF_dev <- reduceDimension(gliomaF_dev, max_components = 2,method = 'DDRTree')
#Trajectory step 3: order cells along the trajectory  
gliomaF_dev <- orderCells(gliomaF_dev)
save(gliomaF_dev,file="data/gliomaF_dev.RData")
pdf("results/pca16res0.3/pseudotime/plot_cell_trajectory_gliomaF.pdf", width = 24, height = 18)
plot_cell_trajectory(gliomaF_dev, color_by = "seurat_clusters")
plot_cell_trajectory(gliomaF_dev, color_by = "Pseudotime")
plot_cell_trajectory(gliomaF_dev,color_by="cell_type")
plot_cell_trajectory(gliomaF_dev,color_by="label")+facet_wrap(~label, nrow = 1)
dev.off()

pdf("results/pca16res0.3/pseudotime/genes.discelltype.gliomaF.pdf", width = 27, height = 18)
to_be_tested <- row.names(subset(fData(monocle_gliomaF_cds),
              gene_short_name %in% c("PTEN", "SOX2", "OLIG1", "MEG3")))
cds_subset_gliomaF <- monocle_gliomaF_cds[to_be_tested,]
plot_genes_jitter(cds_subset_gliomaF,
                  grouping = "cell_type",
                  color_by = "cell_type",
                  nrow= 2,
                  ncol = NULL,
                  plot_trend = TRUE)
dev.off()
pdf("results/pca16res0.3/pseudotime/genes.changeasafunctionofpsedotime.gliomaF.pdf", width = 27, height = 18)
diff_test_res <- differentialGeneTest(gliomaF_dev,fullModelFormulaStr = "~sm.ns(Pseudotime)")
save(diff_test_res,file="data/diff_test_res_gliomaF.Rdata")
diff_test_res[,c("gene_short_name", "pval", "qval")]
diff_test_res_filterbyExp <- subset(diff_test_res, num_cells_expressed >= 10)
#diff_test_res_filterbySig <- subset(diff_test_res_filterbyExp, qval < 10e-10)
sig_gene_names <- row.names(subset(diff_test_res_filterbyExp, qval < 10e-10))
pdf("results/pca16res0.3/pseudotime/heatmap.gliomaF.pdf")
plot_pseudotime_heatmap(gliomaF_dev[sig_gene_names,],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)
dev.off()

pdf("../2single_cell_rnaseq/macrophage/macrophage.pdf", width = 10, height = 10)
DimPlot(glioma.combined.pca16, reduction = "umap",group.by ="celltype",cols=c('Macrophage'='red','Cancer stem cell'='grey','Oligodendrocyte'='grey','Microglia'='grey','Astrocyte'='grey','T cell'='grey','Multilymphoid progenitor cell'='grey','Neural progenitor cell'='grey','Mural cell'='grey','Oligodendrocyte progenitor cell'='grey','Endotheliocyte'='grey'))+ NoLegend()
dev.off()

macro.num<-matrix(data=c(11547, 5952, 0.2401, 0.1238), nrow = 2, ncol = 2, byrow = FALSE, dimnames = NULL)
rownames(macro.num)<-c("primary tumor","recurrent tumor")
colnames(macro.num)<-c("the number of cells","the percentage of cells")
par(mar = c(5, 5, 3, 4)+0.1) #似乎是设置图片位置
bar<-barplot(macro.num[1:2,1],ylim=c(0,12000),ylab="the number of cells",col="blue",col.axis="blue",col.lab="blue",width = 0.1)
par(new=T)
plot(bar,macro.num[1:2,2],axes=F,ylim=c(0,0.3),xlab="",ylab="",col="red",type="o")
#ylim设置不好的话就会看不见折线
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("the percentage of cells",side=4,line=3,col="red")
par(new=T)
legend("top",c('Macrophage'),col=3,lty=1)

csc.num<-matrix(data=c(1042, 4906, 0.0217, 0.1020), nrow = 2, ncol = 2, byrow = FALSE, dimnames = NULL)
rownames(csc.num)<-c("primary tumor","recurrent tumor")
colnames(csc.num)<-c("the number of cells","the percentage of cells")
par(mar = c(5, 5, 3, 4)+0.1) #似乎是设置图片位置
bar<-barplot(csc.num[1:2,1],ylim=c(0,5000),ylab="the number of cells",col="blue",col.axis="blue",col.lab="blue",width = 0.1)
par(new=T)
plot(bar,csc.num[1:2,2],axes=F,ylim=c(0,0.15),xlab="",ylab="",col="red",type="o")
#ylim设置不好的话就会看不见折线
axis(4,col="red",col.ticks="red",col.axis="red")
mtext("the percentage of cells",side=4,line=3,col="red")
par(new=T)
legend("top",c('Cancer stem cell'),col=3,lty=1)


y<-data.frame(date=c('primary tumor','recurrent tumor'),
              primarytumor = c(11547, 5952),
              recurrenttumor=percent(c(24.01%, 12.38%)))
              
              
barnline<-list(
  title = list(text = 'macrophage'),
  tooltip = list(),
  legend = list(data=c('the number of cells','the percentage of cells')),
  xAxis= list(
    type= 'category',
    data= y$date,
    splitLine = list(show=FALSE) # 删掉竖线
  ),
  yAxis= list(
    list(
      type= 'value',
      name= 'the number of cells',
      min= 0,
      max= ceiling(max(y$waterfall/10))*10,
      interval= ceiling(max(y$waterfall/10))*10/5,
      splitLine=list(
        show=FALSE # 删掉横线
      ) 
      #axisLabel= list (formatter= '{value} ml')
    ),
    list(
      type= 'value',
      name= 'the percentage of cells',
      min= 0,
      splitLine=list(
        show=FALSE
      ) ,
      max= 1,
      interval= 0.2 ,
      axisLabel= list(formatter= JS("function(value){return value * 100 + '%';}"))
    )
  ),
  series = list(
    list(
      name="the number of cells",
      type='bar',
      data=y$waterfall,
      itemStyle = list(
        normal=list(color = '#227487', label = list(show = TRUE)) # set bar color and label
      )
    ),
    list(
      name='the percentage of cells',
      type='line',
      yAxisIndex= 1,
      data=y$humid
![excel.jpg](http://upload-images.jianshu.io/upload_images/4006139-cde6c9d7d654bb1c.jpg?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
ity,
      itemStyle = list(
        normal=list( label = list(show = TRUE, formatter= JS("function(c){return Math.floor(c.value  * 10000)/100 + '%';}"))) # label保留小数点后两位
      )
    )
  )
)
echart(barnline, width = 700, height =400)



pdf("../2single_cell_rnaseq/cancerstemcell/cancerstemcell.pdf", width = 10, height = 10)
DimPlot(glioma.combined.pca16, reduction = "umap",group.by ="celltype",cols=c('Macrophage'='grey','Cancer stem cell'='red','Oligodendrocyte'='grey','Microglia'='grey','Astrocyte'='grey','T cell'='grey','Multilymphoid progenitor cell'='grey','Neural progenitor cell'='grey','Mural cell'='grey','Oligodendrocyte progenitor cell'='grey','Endotheliocyte'='grey'))+ NoLegend()
dev.off()

Macrophage"                      "Macrophage"               
 [3] "Oligodendrocyte"                 "Microglia"                      
 [5] "Astrocyte"                       "T cell"                         
 [7] "Multilymphoid progenitor cell"   "Neural progenitor cell"         
 [9] "Mural cell"                      "Oligodendrocyte progenitor cell"
[11] "Endotheliocyte"

macrophage<-subset(x=glioma.combined.pca16,idents="Macrophage")
macrophage_gliomaP<-subset(x=macrophage,subset = label == "gliomaP")
macrophage_gliomaF<-subset(x=macrophage,subset = label == "gliomaF")
cancerstemcell<-subset(x=glioma.combined.pca16,idents="Cancer stem cell")
cancerstemcell_gliomaP<-subset(x=cancerstemcell,subset = label == "gliomaP")
cancerstemcell_gliomaF<-subset(x=cancerstemcell,subset = label == "gliomaF")


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
macrophage.cells <- subset(glioma.combined.pca16, idents = "Macrophage")
Idents(macrophage.cells) <- "label"
avg.macrophage.cells <- as.data.frame(log1p(AverageExpression(macrophage.cells, verbose = FALSE)$RNA))
avg.macrophage.cells$gene <- rownames(avg.macrophage.cells)
write.table(avg.macrophage.cells,"macrophage/avg.macrophage.cells.csv",sep=",",quote=F)
glioma.combined.pca16$celltype.label <- paste(Idents(glioma.combined.pca16), glioma.combined.pca16$label, sep = "_")
glioma.combined.pca16$celltype <- Idents(glioma.combined.pca16)
Idents(glioma.combined.pca16) <- "celltype.label"
immune.response <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP", verbose = FALSE)
head(immune.response, n = 15)
write.table(immune.response,"results/pca16res0.3/degs/macrophage.degs.csv",sep=",")
#CXCL10 M1的marker; S100B在复发中显著下调（促进glioma细胞增殖）；CXCL3在复发中下调（CXCL3是促进炎症的免疫反应的）
pdf("results/pca16res0.3/degs/macrophage.figure1.pdf",width = 20, height = 18)
FeaturePlot(glioma.combined.pca16, features = c("S100B","CXCL3"), split.by = "label", max.cutoff = 3, 
    cols = c("grey", "red"))
dev.off()
pdf("results/pca16res0.3/degs/macrophage.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("S100B","CXCL3"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("results/pca16res0.3/degs/macrophage.figure3.pdf")
genes.to.label = c("S100B","CXCL3")
p1 <- ggplot(avg.macrophage.cells, aes(gliomaF, gliomaP)) + geom_point() + ggtitle("Macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()

pdf("macrophage/TF.SPI1.CEBPB.pdf",width = 20, height = 18)
FeaturePlot(glioma.combined.pca16, features = c("SPI1","CEBPB"), split.by = "label", max.cutoff = 3, 
    cols = c("grey", "red"))
dev.off()
pdf("macrophage/TF.SPI1.CEBPB.fig2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("SPI1","CEBPB"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf("macrophage/TF.CEBPA.NR4A1.pdf",width = 20, height = 18)
FeaturePlot(glioma.combined.pca16, features = c("CEBPA","NR4A1"), split.by = "label", max.cutoff = 3, 
    cols = c("grey", "red"))
dev.off()
pdf("macrophage/TF.CEBPA.NR4A1.fig2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("CEBPA","NR4A1"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf("results/pca16res0.3/degs/macrophage.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("S100B","CXCL3"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("results/pca16res0.3/degs/macrophage.figure3.pdf")
genes.to.label = c("S100B","CXCL3")
p1 <- ggplot(avg.macrophage.cells, aes(gliomaF, gliomaP)) + geom_point() + ggtitle("Macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()





cancerstemcell.cells <- subset(glioma.combined.pca16, idents = "Cancer stem cell")
Idents(cancerstemcell.cells) <- "label"
avg.cancerstemcell.cells <- as.data.frame(log1p(AverageExpression(cancerstemcell.cells, verbose = FALSE)$RNA))
avg.cancerstemcell.cells$gene <- rownames(avg.cancerstemcell.cells)
write.table(avg.macrophage.cells,"cancerstemcell/avg.macrophage.cells.csv",sep=",",quote=F)
glioma.combined.pca16$celltype.label <- paste(Idents(glioma.combined.pca16), glioma.combined.pca16$label, sep = "_")
glioma.combined.pca16$celltype <- Idents(glioma.combined.pca16)
Idents(glioma.combined.pca16) <- "celltype.label"
immune.response <- FindMarkers(glioma.combined.pca16, ident.1 = "Cancer stem cell_gliomaF", ident.2 = "Cancer stem cell_gliomaP", verbose = FALSE)
head(immune.response, n = 15)
write.table(immune.response,"results/pca16res0.3/degs/cancerstemcell.degs.csv",sep=",")
#CXCL10 M1的marker; S100B在复发中显著下调（促进glioma细胞增殖）；CXCL3在复发中下调（CXCL3是促进炎症的免疫反应的）
pdf("results/pca16res0.3/degs/cancerstemcell.figure1.pdf",width = 20, height = 18)
FeaturePlot(glioma.combined.pca16, features = c("VEGFA","SOX4"), split.by = "label", max.cutoff = 3, 
    cols = c("grey", "red"))
dev.off()
pdf("results/pca16res0.3/degs/cancerstemcell.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("VEGFA","SOX4"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("results/pca16res0.3/degs/cancerstemcell.figure3.pdf")
genes.to.label = c("VEGFA","SOX4")
p1 <- ggplot(avg.cancerstemcell.cells, aes(gliomaF, gliomaP)) + geom_point() + ggtitle("cancerstemcell")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()


pdf("cancerstemcell/TF.TEAD4.SOX4.pdf",width = 20, height = 18)
FeaturePlot(glioma.combined.pca16, features = c("TEAD4","SOX4"), split.by = "label", max.cutoff = 3, 
    cols = c("grey", "red"))
dev.off()
pdf("cancerstemcell/TF.TEAD4.SOX4.fig2.pdf",width = 10, height = 18)
plots <- VlnPlot(glioma.combined.pca16, features = c("TEAD4","SOX4"), split.by = "label", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()

pdf("MYB.pdf")
plots <- VlnPlot(glioma.combined.pca16, features = "MYB", split.by = "label", group.by = "celltype.label", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("TMED10.pdf")
plots <- VlnPlot(glioma.combined.pca16, features = "TMED10", split.by = "label", group.by = "celltype.label", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("LARP7.pdf")
plots <- VlnPlot(glioma.combined.pca16, features = "LARP7", split.by = "label", group.by = "celltype.label", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("COMT.pdf")
plots <- VlnPlot(glioma.combined.pca16, features = "COMT", split.by = "label", group.by = "celltype.label", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
MYB  TMED10  LARP7  COMT 




pdf("MYB.distribution.pdf",width=10,height=10)
FeaturePlot(glioma.combined.pca16, features = "MYB")
dev.off()

pdf("TMED10.distribution.pdf",width=10,height=10)
FeaturePlot(glioma.combined.pca16, features = "TMED10")
dev.off()


pdf("LARP7.distribution.pdf",width=10,height=10)
FeaturePlot(glioma.combined.pca16, features = "LARP7")
dev.off()


pdf("COMT.distribution.pdf",width=10,height=10)
FeaturePlot(glioma.combined.pca16, features = "COMT")
dev.off()




require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(macrophage_gliomaP))
names(geneList)<-row.names(macrophage_gliomaP)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(glioma.pca16.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(glioma.pca16.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="data/gene_erichment_results.RData")




glioma.combined.pca16@meta.data$cell.type <- Idents(glioma.combined.pca16)
test <- as.SingleCellExperiment(glioma.combined.pca16)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
immu.se <- DatabaseImmuneCellExpressionData()
sceESC<-LaMannoBrainData('human-es')
sceEmidBrain<-LaMannoBrainData('human-embryo')
sceIPSC<-LaMannoBrainData('human-ips')
#sceESC<-sceESC[,!is.na(sceESC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceESC <- logNormCounts(sceESC)
sceEmidBrain<-sceEmidBrain[,!is.na(sceEmidBrain$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceEmidBrain <- logNormCounts(sceEmidBrain)
sceIPSC<-sceIPSC[,!is.na(sceIPSC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceIPSC <- logNormCounts(sceIPSC)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se, ESC=sceESC, midbrain=sceEmidBrain, iPSC=sceIPSC),
                labels = list(hpca.se$label.main , bpe.se$label.main, sceESC$Cell_type, sceEmidBrain$Cell_type, sceIPSC$Cell_type),
                method = "cluster",
                cluster = test$cell.type)            
table(Anno$pruned.labels)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
#将细胞注释信息重新添加到Seurat对象中去
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(glioma.combined.pca16)
glioma.combined.pca16 <- RenameIdents(glioma.combined.pca16, new.cluster.ids)
head(Idents(glioma.combined.pca16), 5)
levels(Idents(glioma.combined.pca16))
# Visualization
pdf(paste0("results/predicted/figure1.pdf"), width = 24, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label",pt.size = 2)
p2 <- DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE,pt.size = 2)
plot_grid(p1, p2)
dev.off()
pdf(paste0("results/predicted/figure2.pdf"), width = 24, height = 18)
DimPlot(glioma.combined.pca16, reduction = "umap", split.by = "label")
dev.off()
glioma.pca16.markers <- FindAllMarkers(object = glioma.combined.pca16, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(glioma.pca16.markers,"results/predicted/glioma.combined.pca16.markers.csv",sep=",",quote=F)
top30<-glioma.pca16.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_logFC)
write.table(top30,"results/predicted/top30.markers.csv",sep=",",quote=F)





# Run non-linear dimensional reduction (UMAP/tSNE)
glioma.combined.pca16 <- RunUMAP(glioma.combined.pca16, reduction = "pca", dims = 1:18)
# Visualization
pdf(paste0("figures/umap.pca16.res",r,".pdf"), width = 24, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label")
p2 <- DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
pdf(paste0("figures/umap.pca16.res",r,".cluster.pdf"), width = 24, height = 18)
DimPlot(glioma.combined.pca16, reduction = "umap", split.by = "label")
dev.off()
save(glioma.combined.pca16, file="data/autoAnnoglioma.combined.pca16.RData")





prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label))
allsampleprop <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16), glioma.combined.pca16$label)))
write.csv(x = allsampleprop,file = 'anno.allsample.prop.csv',quote = T,row.names = T)
table(Idents(glioma.combined.pca16))
table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
pro <- table(Idents(glioma.combined.pca16),glioma.combined.pca16$label)
write.csv(x =pro,file = 'anno.pro.csv',quote = T,row.names = T)










#Macrophage
pdf("results/pca16res0.3/markergenes/markergenes.Macrophage.pdf",width=15,height=20)
VlnPlot(glioma.combined.pca16, features = c("SPP1","GPNMB","APOC1","FTL","CTSB","S100A9","S100A8","HLA-DQB1","FOLR2","MS4A4A","HLA-DQA1","STAB1","SERPINA1"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Cancer stem cell
pdf("results/pca16res0.3/markergenes/markergenes.CancerStemCell.pdf",width=15,height=20)
VlnPlot(glioma.combined.pca16, features = c("SEC61G","PTN","EGFR","PTPRZ1","C1orf61","SOX4","MEG3","NOVA1","SOX2","GRIA2","GPM6A","CCND2","BEX1"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Oligodendrocyte
pdf("results/pca16res0.3/markergenes/markergenes.Oligodendrocyte.pdf",width=15,height=20)
VlnPlot(glioma.combined.pca16, features = c("PLP1","TF","PPP1R14A","PTGDS","MBP","NKX6-2","MAG","CLDND1","CRYAB","APLP1","ERMN","CLDN11","CNDP1"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Microglia
pdf("results/pca16res0.3/markergenes/markergenes.Microglia.pdf",width=15,height=15)
VlnPlot(glioma.combined.pca16, features = c("HLA-DRB5","CCL3","CCL4","IL1B","CCL3L1","CD83","CCL3L3","TNF","EGR2","NR4A1","CCL4L2"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Astrocyte
pdf("results/pca16res0.3/markergenes/markergenes.Astrocyte.pdf",width=15,height=15)
VlnPlot(glioma.combined.pca16, features = c("CHI3L1","MT3","MT2A","GFAP","RAMP1","IGFBP7","RARRES2","MGST1","SPARCL1"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#T cell
pdf("results/pca16res0.3/markergenes/markergenes.Tcell.pdf",width=15,height=15)
VlnPlot(glioma.combined.pca16, features = c("NKG7","CCL5","IL32","GZMA","CD2","CD52","GZMH","CD3D","GZMB","PRF1"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Multilymphoid progenitor cell
pdf("results/pca16res0.3/markergenes/markergenes.Multilymphoidprogenitorcell.pdf",width=15,height=15)
VlnPlot(glioma.combined.pca16, features = c("HSPB1","HSPA6","PDK4","ZFAND2A","HSPH1","KMT2E-AS1","CDKN1A","FAM46A","GADD45G","OTUD1","DDIT3"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Neural progenitor cell
pdf("results/pca16res0.3/markergenes/markergenes.Neural progenitor cell.pdf",width=15,height=20)
VlnPlot(glioma.combined.pca16, features = c("HIST1H4C","CENPF","TOP2A","MKI67","TUBB","UBE2C","BIRC5","NUSAP1","ASPM","PTTG1","TYMS","KIAA0101","TPX2"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()


#Mural cell
pdf("results/pca16res0.3/markergenes/markergenes.Muralcell.pdf",width=15,height=20)
VlnPlot(glioma.combined.pca16, features = c("DCN","COL3A1","COL1A1","APOD","COL1A2","LUM","MGP","NDUFA4L2","COL4A1","COL4A2","NR2F2","LAMB1","CFH","NID1","PDGFRB"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Oligodendrocyte progenitor cell
pdf("results/pca16res0.3/markergenes/markergenes.Oligodendrocyteprogenitorcell.pdf",width=15,height=10)
VlnPlot(glioma.combined.pca16, features = c("BCAS1","FXYD6","SCRG1","FYN","CDKN1C","SIRT2","RP11-161M6.2","TNR"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()

#Endotheliocyte
pdf("results/pca16res0.3/markergenes/markergenes.Endotheliocyte.pdf",width=15,height=10)
VlnPlot(glioma.combined.pca16, features = c("CLDN5","ITM2A","ANGPT2","VWF","SDPR","EDN1","SLC38A5","GPR116"), group.by = 'celltype', pt.size = 0,ncol=4) 
dev.off()


load("data/glioma.combined.pca16.res0.3.afteranno.RData")
load("data/glioma.pca16.markers.afteranno.RData")
glioma.combined.pca16_gliomaP<-subset(x=glioma.combined.pca16,subset = label == "gliomaP")
glioma.combined.pca16_gliomaF<-subset(x=glioma.combined.pca16,subset = label == "gliomaF")
macrophage.cells <- subset(glioma.combined.pca16, idents = "Macrophage")
Idents(macrophage.cells) <- "label"
avg.macrophage.cells <- as.data.frame(log1p(AverageExpression(macrophage.cells, verbose = FALSE)$RNA))
avg.macrophage.cells$gene <- rownames(avg.macrophage.cells)
glioma.combined.pca16$celltype.label <- paste(Idents(glioma.combined.pca16), glioma.combined.pca16$label, sep = "_")
glioma.combined.pca16$celltype <- Idents(glioma.combined.pca16)
Idents(glioma.combined.pca16) <- "celltype.label"
#immune.response <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP", verbose = FALSE)

## GSVA: https://www.jianshu.com/p/a1dc7380287d?utm_campaign=hugo
#https://www.jianshu.com/p/a1dc7380287d
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library("GSVA")
#remotes::install_version(package = 'Seurat', version = package_version('3.1.4'))
library("Seurat")
packageVersion("Seurat")
# 4.0.0
#DD1 <-  SubsetData(glioma.combined.pca16,
#                   # 提取数据根据的组名,提取两个亚群分群
#                   subset.name = "celltype.label", 
#                   # 提取的组别
#                   accept.value = c("Macrophage_gliomaF","Macrophage_gliomaP"))                 
DD1 <-  subset(glioma.combined.pca16, subset = celltype.label == "Macrophage_gliomaF" | celltype.label == "Macrophage_gliomaP")
#DD1 <-  subset(glioma.combined.pca16, celltype.label == c("Macrophage_gliomaF","Macrophage_gliomaP"))
expr <- as.data.frame(DD1@assays$RNA@data)
expr$ID <- as.matrix(rownames(expr))
s2e <- bitr(expr$ID, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类 转换成ENTREZID
expr <- inner_join(expr,s2e,by=c("ID"="SYMBOL"))
rownames(expr) <- expr$ENTREZID
expr <- expr[,-17501]
expr <- expr[,-17500]
meta <- DD1@meta.data[,c("celltype.label")]
kegggmt <- read.gmt("data/gsea/h.all.v7.2.entrez.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
expr=as.matrix(expr)
kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
save(kegg2,file="results/pca16res0.3/GSVA/macrophage.kegg2.RData")
library(limma)
exprSet <- kegg2
### 差异分析
## 创建分组
group <- factor(meta,levels = c("Macrophage_gliomaP","Macrophage_gliomaF"),ordered = F)
## 分组变成向量，并且限定leves的顺序
## levels里面，把对照组放在前面
## 1.构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)
##2.线性模型拟合
fit <- lmFit(exprSet,design)
##3.贝叶斯检验
fit2 <- eBayes(fit)
#4.输出差异分析结果,其中coef的数目不能操过design的列数
# 此处的2代表的是第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
library('EnhancedVolcano')
adjPvalueCutoff <- 0.05
logFCcutoff <- log2(1.5)
allGeneSets <- topTable(fit2, coef=2, number=Inf)
DEgeneSets <- topTable(fit2, coef=2, number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit2, p.value=adjPvalueCutoff)
summary(res)
pdf("VP.geneset.macrophage.pdf")
DEgeneSets$significant="no"
#DEgeneSets$significant=ifelse(DEgeneSets$logFC>0|DEgeneSets$logFC<0,"up","down")
DEgeneSets$significant=ifelse(DEgeneSets$logFC>0,"up","down")
summary(DEgeneSets$significant)
#ggplot(DEgeneSets,aes(logFC,-1*log10(adj.P.Val)))+geom_point(aes(color =significant)) + xlim(-4,4) + ylim(0,30)+labs(title="Volcanoplot",x="log[2](FC)", y="-log[10](FDR)")+scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)
ggplot(DEgeneSets,aes(logFC,-1*log10(adj.P.Val+0.001)))+geom_point(aes(color =significant)) + xlim(-1,1) + ylim(0,3)+labs(title="Volcanoplot",x="log[2](FC)", y="-log[10](FDR)")+scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+geom_hline(yintercept=1.30103,linetype=4)+geom_vline(xintercept=c(-0.5849625,0.5849625),linetype=4)
dev.off()
#Thus, there are  MSigDB C2 curated pathways that are differentially activated between Macrophage_gliomaP and Macrophage_gliomaF at 0.05 FDR.
allGenes <- topTable(fit2, coef=2, number=Inf)
DEgenes <- topTable(fit2, coef=2, number=Inf,
                    p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit2, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)
pdf("VP.gene.macrophage.pdf")
DEgenes$significant="no"
DEgenes$significant=ifelse(DEgeneSets$logFC>0,"up","down")
summary(DEgenes$significant)
ggplot(DEgenes,aes(logFC,-1*log10(adj.P.Val+0.001)))+geom_point(aes(color =significant)) + xlim(-1,1) + ylim(0,3)+labs(title="Volcanoplot",x="log[2](FC)", y="-log[10](FDR)")+scale_color_manual(values =c("#00ba38","#619cff","#f8766d"))+geom_hline(yintercept=1.30103,linetype=4)+geom_vline(xintercept=c(-0.5849625,0.5849625),linetype=4)
dev.off()
save(allDiff,file="results/pca16res0.3/GSVA/macrophage.gsva.hallmark.results.RData")
write.table(allDiff,"results/pca16res0.3/GSVA/macrophage.gsva.hallmark.results.csv",sep=",")
#根据allDiff找出你自己感兴趣的通路
#up <- c("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_HYPOXIA","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_ANGIOGENESIS","HALLMARK_GLYCOLYSIS","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_MYC_TARGETS_V1","HALLMARK_ADIPOGENESIS","HALLMARK_PROTEIN_SECRETION","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_DNA_REPAIR","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_COAGULATION","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_NOTCH_SIGNALING","HALLMARK_COMPLEMENT","HALLMARK_UV_RESPONSE_UP","HALLMARK_PEROXISOME","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_APOPTOSIS","HALLMARK_HEME_METABOLISM","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_P53_PATHWAY","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_ANDROGEN_RESPONSE","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_HEDGEHOG_SIGNALING","HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_MYC_TARGETS_V2","HALLMARK_APICAL_JUNCTION","HALLMARK_MITOTIC_SPINDLE","HALLMARK_UV_RESPONSE_DN","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_MYOGENESIS","HALLMARK_BILE_ACID_METABOLISM","HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_PANCREAS_BETA_CELLS","HALLMARK_APICAL_SURFACE")
#down <- c("HALLMARK_KRAS_SIGNALING_DN","HALLMARK_SPERMATOGENESIS","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_INTERFERON_GAMMA_RESPONSE")
up <- c("HALLMARK_ANGIOGENESIS","HALLMARK_GLYCOLYSIS","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_ADIPOGENESIS","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
TEST <- up
p <- allDiff
p$ID <- rownames(p) 
q <- p[TEST,]
group1 <- c(rep("Macrophage_gliomaF",45),rep("Macrophage_gliomaP",5)) 
df <- data.frame(ID = q$ID, score = q$t,group=group1 )
# 按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列
head(sortdf)
pdf("results/pca16res0.3/GSVA/macrophage.gsva.figure1.pdf",width = 15, height = 10)
t <- ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity') + 
                    coord_flip() + 
                    theme_bw() + #去除背景色
                    theme(panel.grid =element_blank())+
                    theme(panel.border = element_rect(size = 0.6)) 
t
dev.off()
kegg3 <- as.data.frame(kegg2) 
kegg4 <- kegg3[TEST,]
meta <- DD1@meta.data[,c("celltype.label")]
group <- meta
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(kegg2)

mydata<-kegg4


library(pheatmap)
pdf("results/pca16res0.3/GSVA/macrophage.gsva.figure2.pdf",width = 140,height = 7)
kk2 <- pheatmap(kegg4,
                cluster_rows = F,
                cluster_cols = F,
                annotation_col =annotation_col,
                annotation_legend=T, 
                show_rownames = T,
                show_colnames = F,
                color =colorRampPalette(c("blue", "white","red"))(100),
                cellwidth = 0.5, cellheight = 13,
                fontsize = 10)
print(kk2)
dev.off()

macrophage.allgene <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP")
dim(macrophage.allgene)
write.table(macrophage.de.markers,"macrophage/DEGs.gliomaFandgliomaP.nofilter.csv",sep=",")
macrophage.de.markers <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP", min.pct = 0.5,logfc.threshold = log(1.5),  max.cells.per.ident = 200)
dim(macrophage.de.markers)
head(macrophage.de.markers)
write.table(macrophage.de.markers,"macrophage/DEGs.gliomaFandgliomaP.filter.csv",sep=",")

csc.allgene <- FindMarkers(glioma.combined.pca16, ident.1 = "Cancer stem cell_gliomaF", ident.2 = "Cancer stem cell_gliomaP")
dim(csc.allgene)
write.table(csc.de.markers,"cancerstemcell/DEGs.gliomaFandgliomaP.nofilter.csv",sep=",")
csc.de.markers <- FindMarkers(glioma.combined.pca16, ident.1 = "Cancer stem cell_gliomaF", ident.2 = "Cancer stem cell_gliomaP", min.pct = 0.5,logfc.threshold = log(1.5),  max.cells.per.ident = 200)
dim(csc.de.markers)
head(csc.de.markers)
write.table(csc.de.markers,"cancerstemcell/DEGs.gliomaFandgliomaP.filter.csv",sep=",")


#The results data frame has the following columns :
#p_val : p_val (unadjusted)
#avg_log2FC : log fold-change of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
#pct.1 : The percentage of cells where the feature is detected in the first group
#pct.2 : The percentage of cells where the feature is detected in the second group
#p_val_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.

volcano_names <- ifelse(abs(palmieri_fit_CD$coefficients)>=1,
                        palmieri_fit_CD$genes$SYMBOL, NA)

limma::volcanoplot(palmieri_fit_CD, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)



library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library("GSVA")
DD1 <-  subset(glioma.combined.pca16, subset = celltype.label == "Cancer stem cell_gliomaF" | celltype.label == "Cancer stem cell_gliomaP")
#DD1 <-  SubsetData(glioma.combined.pca16,
#                   # 提取数据根据的组名,提取两个亚群分群
#                   subset.name = "celltype.label", 
#                   # 提取的组别
#                   accept.value = c("Cancer stem cell_gliomaF","Cancer stem cell_gliomaP")) 
expr <- as.data.frame(DD1@assays$RNA@data)
expr$ID <- as.matrix(rownames(expr))
s2e <- bitr(expr$ID, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类 转换成ENTREZID
expr <- inner_join(expr,s2e,by=c("ID"="SYMBOL"))
rownames(expr) <- expr$ENTREZID
expr <- expr[,-5950]
expr <- expr[,-5949]
meta <- DD1@meta.data[,c("celltype.label")]
kegggmt <- read.gmt("data/gsea/h.all.v7.2.entrez.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
expr=as.matrix(expr)
kegg2 <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
save(kegg2,file="results/pca16res0.3/GSVA/cancerstemcell.kegg2.RData")
library(limma)
exprSet <- kegg2
### 差异分析
## 创建分组
group <- factor(meta,levels = c("Cancer stem cell_gliomaP","Cancer stem cell_gliomaF"),ordered = F)
## 分组变成向量，并且限定leves的顺序
## levels里面，把对照组放在前面
## 1.构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)
##2.线性模型拟合
fit <- lmFit(exprSet,design)
##3.贝叶斯检验
fit2 <- eBayes(fit)
#4.输出差异分析结果,其中coef的数目不能操过design的列数
# 此处的2代表的是第二列和第一列的比较
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
save(allDiff,file="results/pca16res0.3/GSVA/cancerstemcell.gsva.results.RData")
write.table(allDiff,"results/pca16res0.3/GSVA/cancerstemcell.gsva.results.csv",sep=",")
up <- c("HALLMARK_MYC_TARGETS_V1","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_DNA_REPAIR","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_MTORC1_SIGNALING","HALLMARK_ADIPOGENESIS","HALLMARK_PROTEIN_SECRETION","HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_NOTCH_SIGNALING","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_GLYCOLYSIS","HALLMARK_E2F_TARGETS","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_MYC_TARGETS_V2","HALLMARK_HYPOXIA","HALLMARK_UV_RESPONSE_UP","HALLMARK_ANDROGEN_RESPONSE","HALLMARK_PEROXISOME","HALLMARK_G2M_CHECKPOINT","HALLMARK_APOPTOSIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_P53_PATHWAY","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HEME_METABOLISM","HALLMARK_MITOTIC_SPINDLE","HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_COAGULATION","HALLMARK_UV_RESPONSE_DN","HALLMARK_COMPLEMENT","HALLMARK_HEDGEHOG_SIGNALING","HALLMARK_APICAL_JUNCTION","HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_ANGIOGENESIS","HALLMARK_PANCREAS_BETA_CELLS","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_BILE_ACID_METABOLISM","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE")
down <- c("HALLMARK_MYOGENESIS","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_APICAL_SURFACE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_ESTROGEN_RESPONSE_EARLY","HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_SPERMATOGENESIS","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_KRAS_SIGNALING_DN")
TEST <- c(up,down)
p <- allDiff
p$ID <- rownames(p) 
q <- p[TEST,]
group1 <- c(rep("Cancer stem cell_gliomaF",41),rep("Cancer stem cell_gliomaP",9)) 
df <- data.frame(ID = q$ID, score = q$t,group=group1 )
# 按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列
head(sortdf)
pdf("results/pca16res0.3/GSVA/cancerstemcell.gsva.figure1.pdf",width = 15, height = 10)
t <- ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity') + 
                    coord_flip() + 
                    theme_bw() + #去除背景色
                    theme(panel.grid =element_blank())+
                    theme(panel.border = element_rect(size = 0.6)) 
t
dev.off()





library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
load("data/glioma.combined.pca16.res0.3.afteranno.RData")
load("data/glioma.pca16.markers.afteranno.RData")
glioma.combined.pca16$celltype.label <- paste(Idents(glioma.combined.pca16), glioma.combined.pca16$label, sep = "_")
glioma.combined.pca16$celltype <- Idents(glioma.combined.pca16)
Idents(glioma.combined.pca16) <- "celltype.label"
macrophage.de.markers.filter <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP", min.pct = 0.5,logfc.threshold = log2(1.5),   base = 2, fc.name = NULL,max.cells.per.ident = 200)
macrophage.de.markers.nofilter <- FindMarkers(glioma.combined.pca16, ident.1 = "Macrophage_gliomaF", ident.2 = "Macrophage_gliomaP", logfc.threshold = 0,min.pct= 0.5,   base = 2, fc.name = NULL,max.cells.per.ident = 200)
dim(macrophage.de.markers.filter)
dim(macrophage.de.markers.nofilter)
data<-subset(macrophage.de.markers.nofilter, macrophage.de.markers.nofilter$p_val_adj <= 0.05)
newdata<-subset(macrophage.de.markers.filter, macrophage.de.markers.filter$p_val_adj <= 0.05)
dim(data)
dim(newdata)
write.table(newdata,"macrophage/DEGs.csv",quote=F,sep=",")
pdf("macrophage/vpplot.pdf")
plot(data$avg_log2FC,-log2(data$p_val_adj))
points(newdata$avg_log2FC,-log2(newdata$p_val_adj),col="red")
text(newdata$avg_log2FC, -log2(newdata$p_val_adj)+5, labels=rownames(newdata), cex=1.0, col='red')
abline(v=c(-log2(1.5),log2(1.5)),lty=2,lwd = 2,col="red")
dev.off()
csc.de.markers.filter <- FindMarkers(glioma.combined.pca16, ident.1 = "Cancer stem cell_gliomaF", ident.2 = "Cancer stem cell_gliomaP", min.pct = 0.5,logfc.threshold = log2(1.5),   base = 2, fc.name = NULL,max.cells.per.ident = 200)
csc.de.markers.nofilter <- FindMarkers(glioma.combined.pca16, ident.1 = "Cancer stem cell_gliomaF", ident.2 = "Cancer stem cell_gliomaP", logfc.threshold = 0,min.pct= 0.5,   base = 2, fc.name = NULL,max.cells.per.ident = 200)
dim(csc.de.markers.filter)
dim(csc.de.markers.nofilter)
data<-subset(csc.de.markers.nofilter, csc.de.markers.nofilter$p_val_adj <= 0.05)
newdata<-subset(csc.de.markers.filter, csc.de.markers.filter$p_val_adj <= 0.05)
dim(data)
dim(newdata)
write.table(newdata,"cancerstemcell/DEGs.csv",quote=F,sep=",")
pdf("cancerstemcell/vpplot.pdf")
plot(data$avg_log2FC,-log2(data$p_val_adj))
points(newdata$avg_log2FC,-log2(newdata$p_val_adj),col="red")
text(newdata$avg_log2FC, -log2(newdata$p_val_adj)+5, labels=rownames(newdata), cex=1.0, col='red')
abline(v=c(-log2(1.5),log2(1.5)),lty=2,lwd = 2,col="red")
dev.off()

















library(Seurat)
library(tidyverse)
dir.create("CellTalk")
#scRNA <- readRDS("scRNA.rds")
#提取HNSCC肿瘤样本HNC01TIL
sp1 <- glioma.combined.pca16[,str_detect(colnames(glioma.combined.pca16),'HNC01TIL')]
sp1_counts <- as.matrix(sp1@assays$RNA@data)
sp1_counts <- data.frame(Gene=rownames(sp1_counts), sp1_counts)
sp1_meta <- data.frame(Cell=rownames(sp1@meta.data), cell_type=sp1@meta.data$celltype_Monaco)
write.table(sp1_counts, "CellTalk/sp1_counts.txt", row.names=F, sep='\t')
write.table(sp1_meta, "CellTalk/sp1_meta.txt", row.names=F, sep='\t')


glioma.combined.pca16_gliomaP<-subset(x=glioma.combined.pca16,subset = label == "gliomaP")
glioma.combined.pca16_gliomaF<-subset(x=glioma.combined.pca16,subset = label == "gliomaF")
write.table(as.matrix(glioma.combined.pca16_gliomaP@assays$RNA@data), 'cellphonedb_gliomaP_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(glioma.combined.pca16_gliomaP@meta.data), glioma.combined.pca16_gliomaP@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_gliomaP_meta.txt', sep='\t', quote=F, row.names=F)

write.table(as.matrix(glioma.combined.pca16_gliomaF@assays$RNA@data), 'cellphonedb_gliomaF_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(glioma.combined.pca16_gliomaF@meta.data), glioma.combined.pca16_gliomaF@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_gliomaF_meta.txt', sep='\t', quote=F, row.names=F)

gliomaP='/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/CellTalk/cellphone_gliomaP/' ##  outs 文件放在这里了。
library(psych)
library(qgraph)
library(igraph)
netf<- "count_network.txt"
mynet <- read.delim(paste0(gliomaP,"count_network.txt"), check.names = FALSE)
head(mynet)
net<- graph_from_data_frame(mynet)
pdf("cellphone_gliomaP/net.fig1.pdf")
plot(net)
dev.off()
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$count/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
     

gliomaF='/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/CellTalk/cellphone_gliomaF/' ##  outs 文件放在这里了。
library(psych)
library(qgraph)
library(igraph)
netf<- "count_network.txt"
mynet <- read.delim(paste0(gliomaF,"count_network.txt"), check.names = FALSE)
head(mynet)
net<- graph_from_data_frame(mynet)
setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/CellTalk/cellphone_gliomaF")
pdf("net.fig1.pdf")
plot(net)
dev.off()
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$count/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
     
     
     
     
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
dir.create("subcluster")
##提取细胞子集
Cells.sub.tme <- subset(glioma.combined.pca16@meta.data, celltype==c("Macrophage","T cell","Multilymphoid progenitor cell"))
scRNAsub.tme <- subset(glioma.combined.pca16, cells=row.names(Cells.sub.tme))
#提重新降维聚类
#因为再聚类的细胞之间差异比较小，所以聚类函数FindClusters()控制分辨率的参数建议调高到resolution = 0.7。
##PCA降维
scRNAsub.tme <- FindVariableFeatures(scRNAsub.tme, selection.method = "vst", nfeatures = 2000)
scale.genes.tme <-  rownames(scRNAsub.tme)
scRNAsub.tme <- ScaleData(scRNAsub.tme, features = scale.genes.tme)
scRNAsub.tme <- RunPCA(scRNAsub.tme, features = VariableFeatures(scRNAsub.tme))
pdf("subcluster/Determine.tme.pcnumber.pdf")
ElbowPlot(scRNAsub.tme, ndims=20, reduction="pca")
dev.off()
pc.num=1:12
##细胞聚类
scRNAsub.tme <- FindNeighbors(scRNAsub.tme, dims = pc.num) 
scRNAsub.tme <- FindClusters(scRNAsub.tme, resolution = 0.7)
table(scRNAsub.tme@meta.data$seurat_clusters)
metadata <- scRNAsub.tme@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'subcluster/tme.cell_cluster.csv',row.names = F)
##非线性降维
#tSNE
scRNAsub.tme = RunTSNE(scRNAsub.tme, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.tme, 'tsne')
write.csv(embed_tsne,'subcluster/tme.embed_tsne.csv')
plot1 = DimPlot(scRNAsub.tme, reduction = "tsne") 
ggsave("subcluster/tme.tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("subcluster/tme.tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub.tme <- RunUMAP(scRNAsub.tme, dims = pc.num)
embed_umap <- Embeddings(scRNAsub.tme, 'umap')
write.csv(embed_umap,'subcluster/tme.embed_umap.csv') 
plot2 = DimPlot(scRNAsub.tme, reduction = "umap") 
ggsave("subcluster/tme.UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("subcluster/tme.UMAP.png", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("subcluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("subcluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
diff.wilcox = FindAllMarkers(scRNAsub.tme)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "subcluster/tme.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "subcluster/tme.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.tme,file="subcluster/scRNAsub.tme.RData")
##细胞类型鉴定
library(SingleR)
refdata <- MonacoImmuneData()
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
immu.se <- DatabaseImmuneCellExpressionData() 
testdata <- GetAssayData(scRNAsub.tme, slot="data")
clusters <- scRNAsub.tme@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, 
										 ref = list(REF=refdata, HP = hpca.se , BP = bpe.se, IM=immu.se),
										 labels = list(refdata$label.fine, hpca.se$label.fine , bpe.se$label.fine, immu.se$label.fine), 
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts",de.method="wilcox")           
table(cellpred$labels)
pdf("subcluster/test.pdf", width=18 ,height=9)
plotScoreHeatmap(cellpred)
dev.off()
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"subcluster/tme.celltype_singleR.csv",row.names = F)
scRNAsub.tme@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
scRNAsub.tme@meta.data[which(scRNAsub.tme@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub.tme, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsub.tme, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("subcluster/tme.tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("subcluster/tme.UMAP_celltype.pdf", p2, width=10 ,height=6)
ggsave("subcluster/tme.celltype.pdf", p3, width=10 ,height=5)
ggsave("subcluster/tme.celltype.png", p3, width=10 ,height=5)

new.cluster.ids<-c("Macrophage:monocyte-derived:IL-4/Dex/cntrl", "Macrophages", "Macrophage:monocyte-derived:IL-4/Dex/TGFb", "Macrophage:monocyte-derived:IL-4/Dex/TGFb", "Macrophage:monocyte-derived:IL-4/Dex/TGFb", "Macrophages M1", "Astrocyte:Embryonic_stem_cell-derived", "CD8+ Tcm", "Macrophages M1", "Astrocyte:Embryonic_stem_cell-derived", "NK_cell", "Monocytes", "Myeloid dendritic cells")
names(new.cluster.ids) <- levels(scRNAsub.tme)
scRNAsub.tme <- RenameIdents(scRNAsub.tme, new.cluster.ids)
scRNAsub.tme <- RunUMAP(scRNAsub.tme, dims = 1:16)
save(scRNAsub.tme,file="data/scRNAsub.tme.afteranno.RData")
scRNAsub.tme.markers <- FindAllMarkers(object = scRNAsub.tme, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.tme.markers,"subcluster/scRNAsub.tme.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.tme.markers,file="subcluster/scRNAsub.tme.markers.afteranno.RData")
pdf("subcluster/scRNAsub.tme.integrate.pdf", width = 36, height = 18)
p1 <- DimPlot(scRNAsub.tme, reduction = "umap", group.by = "label")
p2<-DimPlot(scRNAsub.tme, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
prop.table(table(Idents(scRNAsub.tme), scRNAsub.tme$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.tme), scRNAsub.tme$label)))
prop.table(table(Idents(glioma.combined.pca16)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(glioma.combined.pca16))))
write.csv(x = allsampleprop.each,file = 'results/pca16res0.3/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = 'results/pca16res0.3/anno.allsample.total.prop.csv',quote = T,row.names = T)


pdf("results/pca16res0.3/umap.pca16.res0.3.integrate.pdf", width = 36, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label")
p2<-DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


##细胞周期归类 
merged<- CellCycleScoring(object = glioma.combined.pca16, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes) 
head(x = glioma.combined.pca16@meta.data) 
pdf("cellcycling.pdf")
DimPlot(merged,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()


#conda environment: R40
library(scater)
# install SeuratDisk from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/seurat-disk', ref = 'develop')
library(SeuratDisk)
library(patchwork)
library(Seurat)
library(hdf5r)
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("scater")
library(scater)
library(velocyto.R)
setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities")
load("../data/glioma.combined.pca16.res0.3.afteranno.RData")
DefaultAssay(glioma.combined.pca16)<-'RNA'
glioma.combined.pca16
# seurat对象转换为loom文件
sdata.loom <- as.loom(x = glioma.combined.pca16, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/glioma.combined.pca16.loom", verbose = FALSE)
#sdata.loom <- as.loom(x = glioma.combined.pca16, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/sdata.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)

glioma.combined.pca16_gliomaP<-subset(x=glioma.combined.pca16,subset = label == "gliomaP")
#sdata.loom_gliomaP <- as.loom(x = glioma.combined.pca16_gliomaP, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/sdata_gliomaP.loom", verbose = FALSE)
sdata.loom_gliomaP <- as.loom(x = glioma.combined.pca16_gliomaP, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/gliomaP.loom", verbose = FALSE)
glioma.combined.pca16_gliomaF<-subset(x=glioma.combined.pca16,subset = label == "gliomaF")
#sdata.loom_gliomaF <- as.loom(x = glioma.combined.pca16_gliomaF, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/sdata_gliomaF.loom", verbose = FALSE)
sdata.loom_gliomaF <- as.loom(x = glioma.combined.pca16_gliomaF, filename = "/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/RNAvelocities/gliomaF.loom", verbose = FALSE)

# Always remember to close loom files when done

library(velocyto.R)
ldat <- read.loom.matrices("glioma.combined.pca16.loom")
ldat_gliomaP <- read.loom.matrices("sdata_gliomaP.loom")
ldat_gliomaF <- read.loom.matrices("sdata_gliomaF.loom")


allen_reference <- readRDS("allen_cortex.rds")
testdata <- as.loom(x = allen_reference, filename = "allen_reference.loom", verbose = FALSE)
ldat <- read.loom.matrices("allen_reference.loom")



#导入包
library(velocyto.R)

input_loom <- "gliomaP.loom"
ldat <- read.loom.matrices(input_loom)
ldat <- ReadVelocity(file = input_loom)
#使用剪切位点的表达量作为输入
emat <- ldat$spliced
#做直方图查看数据分步
hist(log10(colSums(emat)),col='wheat',xlab='cell size')


setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/TSCAN")
library(TSCAN)
data(lpsdata)
head(lpsdata)
dim(lpsdata)
class(lpsdata)
load("../data/glioma.combined.pca16.res0.3.afteranno.RData")
DefaultAssay(glioma.combined.pca16)<-'RNA'
glioma.combined.pca16_gliomaP<-subset(x=glioma.combined.pca16,subset = label == "gliomaP")
save(glioma.combined.pca16_gliomaP,file="../data/gliomaP.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaP$celltype<-Idents(glioma.combined.pca16_gliomaP)
glioma.combined.pca16_gliomaF<-subset(x=glioma.combined.pca16,subset = label == "gliomaF")
save(glioma.combined.pca16_gliomaF,file="../data/gliomaF.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaF$celltype<-Idents(glioma.combined.pca16_gliomaF)


seuratdf<-as.matrix(glioma.combined.pca16_gliomaF@assays$RNA@counts)
head(seuratdf)
dim(seuratdf)
rownames(lpsdata)
class(seuratdf)
procdata <- preprocess(seuratdf,cvcutoff = 0)
dim(procdata)
lpsmclust <- exprmclust(procdata)
pdf("PCA_gliomaF.pdf")
plotmclust(lpsmclust,show_cell_names = F)
dev.off()
lpsorder <- TSCANorder(lpsmclust)
lpsorder
diffval <- difftest(procdata,lpsorder)
#Selected differentially expressed genes under qvlue cutoff of 0.05
head(row.names(diffval)[diffval$qval < 0.05])
[1] "RPL22"    "RPL11"    "SH3BGRL3" "CD52"     "LAPTM5"   "YBX1"    
STAT2expr <- log2(seuratdf["WNT6A",]+1)
STAT2expr <- log2(seuratdf["PARK7",]+1)
pdf("testDEGs_gliomaF.pdf")
singlegeneplot(STAT2expr, TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))
dev.off()
TSCANui()# 打开交互界面



setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/Cellrank")
library(sceasy)
library(reticulate)
library(Seurat)
use_condaenv('EnvironmentName')
loompy <- reticulate::import('loompy')
load("../data/gliomaF.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaF$celltype<-Idents(glioma.combined.pca16_gliomaF)
sceasy::convertFormat(glioma.combined.pca16_gliomaF, from="seurat", to="anndata",
                       outFile='glioma.combined.pca16_gliomaF.h5ad')
load("../data/gliomaP.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaP$celltype<-Idents(glioma.combined.pca16_gliomaP)
sceasy::convertFormat(glioma.combined.pca16_gliomaP, from="seurat", to="anndata",
                       outFile='glioma.combined.pca16_gliomaP.h5ad')
import cellrank as cr
adata = cr.datasets.pancreas_preprocessed("glioma.combined.pca16_gliomaF.h5ad")
adata
k = cr.tl.transition_matrix(
    adata, weight_connectivities=0.2, softmax_scale=4, show_progress_bar=False
)
g = cr.tl.estimators.GPCCA(k)
g.fit(n_lineages=3, cluster_key="celltype", compute_absorption_probabilities=True)
g.plot_absorption_probabilities()



rm(list=ls())
options(stringsAsFactors=F)
library(data.table)
library(densityClust)
library(destiny)
library(rgl)
#source('functios.R')
source('colorPalette.R')
#data=females_data
#destinyObj<-as.ExpressionSet(as.data.frame(t(data)))
#destinyObj$condition<-factor(condition)
#dm<-DiffusionMap(destinyObj,sigma,rotate=TRUE)
load("../data/gliomaF.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaF$celltype<-Idents(glioma.combined.pca16_gliomaF)
#functions.R
function(data=data,condition=condition,sigma="local"){
	destinyObj<-as.ExpressionSet(as.data.frame(t(data)))
	destinyObj$condition<-factor(condition)
	dm<-DiffusionMap(destinyObj,sigma,rotate=TRUE)
	return(dm)
}
#compute the diffusion map
female_dm<-run_diffMap(
	female_data,
	female_clustering,
	sigma=15
)

female_dm<-run_diffMap(
	glioma.combined.pca16_gliomaF,
	female_clustering,
	sigma=15
)


save(female_dm,female_clustering,female_stages,file='diffusionMap_ouput.Rdata')

#plot the eigen values per diffusion component (similar to screenplot)
plot_eigenVal(
	dm=female_dm
)



#https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/slingshot")
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)
load("../data/gliomaF.combined.pca16.res0.3.afteranno.RData")
glioma.combined.pca16_gliomaF$celltype<-Idents(glioma.combined.pca16_gliomaF)
write.table(Idents(glioma.combined.pca16_gliomaF),"glioma.combined.pca16_gliomaF.csv",,sep=",",quote=F,col.names=F)


gliomaF.sce <- as.SingleCellExperiment(glioma.combined.pca16_gliomaF)
geneFilter <- apply(assays(gliomaF.sce)$counts,1,function(x){
    sum(x >= 3) >= 10
})
gliomaF.sce <- gliomaF.sce[geneFilter, ]
p1 <- plotExpression(gliomaF.sce, features = "SAMD11", x = "ident") + theme(axis.text.x = element_text(angle = 45, 
    hjust = 1))
p2 <- plotPCA(gliomaF.sce, colour_by = "ident")
pdf("figure1.pdf")
p1 + p2
dev.off()
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(gliomaF.sce)$norm <- FQnorm(assays(gliomaF.sce)$counts)
pca <- prcomp(t(log1p(assays(gliomaF.sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
pdf("figure2.pdf")
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
dev.off()
library(uwot)
rd2 <- umap(t(log1p(assays(gliomaF.sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
pdf("figure3.pdf")
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
dev.off()
reducedDims(gliomaF.sce) <- SimpleList(PCA = rd1, UMAP = rd2)
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(gliomaF.sce)$GMM <- cl1
library(RColorBrewer)
pdf("figure4.pdf")
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
dev.off()
cl2 <- kmeans(rd1, centers = 4)$cluster
colData(gliomaF.sce)$kmeans <- cl2
pdf("figure5.pdf")
plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)
dev.off()
library(slingshot)
sim <- slingshot(gliomaF.sce, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(sim$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
pdf("figure6.pdf")
plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black')
dev.off()
pdf("figure7.pdf")
plot(reducedDims(sim)$PCA, col = brewer.pal(9,'Set1')[sim$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')
dev.off()


library(tradeSeq)
# fit negative binomial GAM
sim <- fitGAM(sim)
# test for dynamic expression
ATres <- associationTest(sim)
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sim$slingPseudotime_1, na.last = NA)
heatdata <- assays(sim)$counts[topgenes, pst.ord]
heatclus <- sim$GMM[pst.ord]
pdf("figure8.pdf")
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
dev.off()

library(slingshot, quietly = FALSE)
data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sim)$counts,1,function(x){
    sum(x >= 3) >= 10
})
sim <- sim[geneFilter, ]
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)



# 表达矩阵
# HVGs
library("BisqueRNA")
library("destiny")
setwd('/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/slingshot')
load("../data/glioma.combined.pca16.res0.3.afteranno.RData")
DefaultAssay(glioma.combined.pca16)="RNA"
glioma.combined.pca16$celltype<-Idents(glioma.combined.pca16)
glioma.combined.pca16 <- FindVariableFeatures(glioma.combined.pca16, selection.method = "vst", nfeatures = 2000)
seurat4_HVGs <- VariableFeatures(glioma.combined.pca16)
glioma.combined.pca16 <- glioma.combined.pca16[rownames(glioma.combined.pca16) %in% as.vector(seurat4_HVGs),]
data<-SeuratToExpressionSet(glioma.combined.pca16,  delimiter='-',  position=2,  version = "v3")
destinyObj <- data
# 2个发育时期获取
head(colnames(glioma.combined.pca16))
stage <- read.csv('glioma.combined.pca16.label.csv')
female_stages=stage$x
stage <- read.csv('glioma.combined.pca16.label.csv')
female_stages=stage$x
names(female_stages) <-glioma.combined.pca16$label
table(female_stages)
# 11个cluster获取
cluster <- read.csv('glioma.combined.pca16.csv')
female_clustering=cluster$x
names(female_clustering)=rownames(cluster)
table(female_clustering)
conditions=female_clustering
destinyObj$condition <- factor(conditions)
sigma=15
dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)
save(dm,"dm.RData")





setwd('/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/slingshot')
load("../data/glioma.combined.pca16.res0.3.afteranno.RData")
DefaultAssay(glioma.combined.pca16)="RNA"
glioma.combined.pca16$celltype<-Idents(glioma.combined.pca16)
glioma.combined.pca16 <- FindVariableFeatures(glioma.combined.pca16, selection.method = "vst", nfeatures = 2000)
seurat4_HVGs <- VariableFeatures(glioma.combined.pca16)
glioma.combined.pca16 <- glioma.combined.pca16[rownames(glioma.combined.pca16) %in% as.vector(seurat4_HVGs),]

write.table(Idents(glioma.combined.pca16),"glioma.combined.pca16.csv",sep=",",quote=F)
females_data<-GetAssayData(object = glioma.combined.pca16, slot = "counts")
write.table(females_data,"glioma.combined.pca16.matrix.csv",sep=",",quote=F)
mymatrix<-read.table("glioma.combined.pca16.matrix.csv",sep=",",header=T,row.names=1)

females_data <- FindVariableFeatures(females, selection.method = "vst", nfeatures = 2000)
dim(females_data)
females_data <- log(females_data+1)
females_data[1:4,1:4]
save(females_data,file = 'females_hvg_matrix.Rdata')

setwd("/home/liuzhe/Wanglab/zhangyu/Glioma_human/2single_cell_rnaseq/slingshot")
load("../data/glioma.combined.pca16.res0.3.afteranno.RData")
DefaultAssay(glioma.combined.pca16)='RNA'
library('matrixStats')
library('statmod')
library("ggplot2")
mymatrix<-read.table("glioma.combined.pca16.matrix.csv",sep=",",header=T,row.names=1)
if(T){
  getMostVarGenes <- function(
    data=data,              # RPKM matrix
    fitThr=1.5,             # Threshold above the fit to select the HGV
    minMeanForFit=1         # Minimum mean gene expression level
  ){
    # data=females;fitThr=2;minMeanForFit=1
    # Remove genes expressed in no cells
    data_no0 <- as.matrix(
      data[rowSums(data)>0,]
    )
    # Compute the mean expression of each genes
    meanGeneExp <- rowMeans(data_no0)
    names(meanGeneExp)<- rownames(data_no0)

    # Compute the squared coefficient of variation
    varGenes <- rowVars(data_no0)
    cv2 <- varGenes / meanGeneExp^2

    # Select the genes which the mean expression is above the expression threshold minMeanForFit
    useForFit <- meanGeneExp >= minMeanForFit

    # Compute the model of the CV2 as a function of the mean expression using GLMGAM
    fit <- glmgam.fit( cbind( a0 = 1,
                              a1tilde = 1/meanGeneExp[useForFit] ),
                       cv2[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])

    # Get the highly variable gene counts and names
    fit_genes <- names(meanGeneExp[useForFit])
    cv2_fit_genes <- cv2[useForFit]
    fitModel <- fit$fitted.values
    names(fitModel) <- fit_genes
    HVGenes <- fitModel[cv2_fit_genes>fitModel*fitThr]
    print(length(HVGenes))

    # Plot the result
    plot_meanGeneExp <- log10(meanGeneExp)
    plot_cv2 <- log10(cv2)
    plotData <-  data.frame(
      x=plot_meanGeneExp[useForFit],
      y=plot_cv2[useForFit],
      fit=log10(fit$fitted.values),
      HVGenes=log10((fit$fitted.values*fitThr))
    )
    p <- ggplot(plotData, aes(x,y)) +
      geom_point(size=0.1) +
      geom_line(aes(y=fit), color="red") +
      geom_line(aes(y=HVGenes), color="blue") +
      theme_bw() +
      labs(x = "Mean expression (log10)", y="CV2 (log10)")+
      ggtitle(paste(length(HVGenes), " selected genes", sep="")) +
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        legend.position= "none",
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        aspect.ratio=1
      )+
      scale_color_manual(
        values=c("#595959","#5a9ca9")
      )
    print(p)

    # Return the RPKM matrix containing only the HVG
    HVG <- data_no0[rownames(data_no0) %in% names(HVGenes),]
    return(HVG)
  }
}
females<-mymatrix
females_data <- getMostVarGenes(females, fitThr=2)

#females_data <- FindVariableFeatures(glioma.combined.pca16, selection.method = "vst", nfeatures = 2000)
data<-SeuratToExpressionSet(glioma.combined.pca16,  delimiter='-',  position=2,  version = "v3")


# 2个发育时期获取
head(colnames(glioma.combined.pca16))
names(female_stages) <-glioma.combined.pca16$label
write.table(glioma.combined.pca16$label,"glioma.combined.pca16.label.csv",,sep=",",quote=F)
stage <- read.csv('glioma.combined.pca16.label.csv')
female_stages=stage$x
names(female_stages)=rownames(stage)
table(female_stages)

# 11个cluster获取
cluster <- read.csv('glioma.combined.pca16.csv')
female_clustering=cluster$x
names(female_clustering)=rownames(cluster)
table(female_clustering)

# 1.2 进行DiffusionMap
female_dm <- run_diffMap(
  females_data, 
  female_clustering,
  sigma=15
)
# 这个包装的函数其实做了下面几行代码的事情
condition=female_clustering
sigma=15
library("BisqueRNA")
library("destiny")
destinyObj <- mymatrix
destinyObj$condition <- factor(condition)
dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)

save(female_dm,females_data,female_clustering,female_stages,
     file = 'diffusionMap_output.Rdata')

data=females_data
condition=female_clustering
sigma=15
destinyObj <- as.ExpressionSet(as.data.frame(t(mymatrix)))
destinyObj$condition <- factor(condition)
dm <- DiffusionMap(destinyObj, sigma, rotate = TRUE)

save(female_dm,females_data,female_clustering,female_stages,
     file = 'diffusionMap_output.Rdata')







# 1.3 作图探索
plot_eigenVal(
  dm=female_dm
)



# 探索4个分群
female_clusterPalette <- c(
  C1="#560047",
  C2="#a53bad", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)

plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_clustering, 
  colour=female_clusterPalette
)
# 探索6个发育时间
female_stagePalette=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_stages, 
  colour=female_stagePalette
)

# 2 进行Slingshot

dm=female_dm
dim=c(1:4)
condition=factor(female_clustering)

data <- data.frame(
  dm@eigenvectors[,dim]
)

female_lineage <- slingshot(
  data, 
  condition, 
  start.clus = "C1", 
  end.clus=c("C2", "C4"),
  maxit=100000,
  shrink.method="cosine"
  # shrink.method="tricube"
)
# 看下结果
> female_lineage
class: SlingshotDataSet 

 Samples Dimensions
     563          4

lineages: 2 
Lineage1: C1  C3  C4  
Lineage2: C1  C2  

curves: 2 
Curve1: Length: 1.3739  Samples: 453.62
Curve2: Length: 0.74646 Samples: 312.73
female_pseudotime <- get_pseudotime(female_lineage, wthres=0.9)
rownames(female_pseudotime) <- colnames(females)