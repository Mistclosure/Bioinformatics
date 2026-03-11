library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(WGCNA)
library(reshape2)
library(stringr)
setwd('/mnt/disk1/qiuzerui/downloads/Rloop')
# 加载C部分保存的恶性细胞子集
load("Malignant.Rdata")
# 2. 核心修改：由于前置步骤使用了 SCTransform，必须从 SCT assay 提取 data 图层
# 这样提取出的 datadf 才带有基因名（行名）和细胞ID（列名），解决 dimnames 错误
pbmc1 <- NormalizeData(pbmc1, assay = "RNA")
datadf <- as.matrix(GetAssayData(pbmc1, assay = "RNA", layer = "data"))
# -----------------------------------------------------

idd1 <- pbmc1@meta.data
# 获取当前恶性细胞子集的聚类信息（来自C部分的分群）
Inter.id1 <- cbind(rownames(idd1), as.character(idd1$seurat_clusters))
rownames(Inter.id1) <- rownames(idd1)
colnames(Inter.id1) <- c("CellID", "Celltype")
Inter.id1 <- as.data.frame(Inter.id1)
head(Inter.id1)

set.seed(9)
# 修正：直接使用 pbmc1 的元数据进行合并
A = merge(pbmc1@meta.data, Inter.id1, by = 0)
Inter.id1 = A[, c(1, ncol(A))] # 动态锁定最后一列作为 Celltype
colnames(Inter.id1) <- c("CellID", "Celltype")
rownames(Inter.id1) = Inter.id1[, 1]

# 此时 datadf 已具备列名，以下切片操作将正常运行
Inter1 <- datadf[, Inter.id1$CellID]
Inter2 <- as.matrix(Inter1)
Inter.id1$Celltype = as.factor(Inter.id1$Celltype)
Inter2[1:4,1:4]

x = Inter.id1$Celltype
pseudocell.size = 30
new_ids_list1 = list()
table(Inter.id1$Celltype)

# 伪细胞构建流程
for (i in 1:length(levels(Inter.id1$Celltype))) {
  cluster_id = levels(Inter.id1$Celltype)[i]
  cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)
  new_ids_list1[[i]] <- pseudo_ids
}

new_ids <- unlist(new_ids_list1)
new_ids <- as.data.frame(new_ids)
head(new_ids)

new_ids_length <- table(new_ids)
new_ids_length
new_colnames <- rownames(new_ids)
colnames(datadf)

all.data <- datadf[, as.character(new_colnames)]
all.data <- t(all.data)
new.data2 <- aggregate(all.data, list(name = new_ids[, 1]), FUN = mean)
rownames(new.data2) <- new.data2$name
new.data2 <- new.data2[, -1]

new_ids_length <- as.matrix(new_ids_length)
new_good_ids <- as.matrix(new_ids_length)
result <- t(new.data2)[, rownames(new_ids_length)]
dim(result)

# --- 读取 Rloop.csv 
rloop = read.csv("Rloop.csv", header = T, check.names = F, row.names = 1)
Cluster1 <- result[intersect(rownames(rloop), rownames(result)), ]
dim(Cluster1)

# --- WGCNA 网络构建 ---
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType == "pearson", cor, bicor)
corFnc
maxPOutliers = ifelse(corType == "pearson", 1, 0.05)
robustY = ifelse(corType=="pearson", T, F)
dataExpr <- as.matrix(t(Cluster1))

gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
gsg$goodSamples
if (!gsg$allOK) { 
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, networkType = "signed", verbose = 5)

par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

power = sft$powerEstimate
if(is.na(power)) power = 12 # 容错处理
softPower = power
softPower

cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, corType = corType,
                       maxPOutliers = maxPOutliers, loadTOMs = TRUE,
                       saveTOMFileBase = paste0("dataExpr", ".tom"),
                       verbose = 3)
table(net$colors)

# 模块颜色与可视化
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleColors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
head(MEs_col)

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(3, 3, 2, 4),
                      marHeatmap = c(3, 4, 2, 2), plotDendrograms = T, xLabelsAngle = 90)

# ====== 恢复原代码缺失的模块表达可视化及TOM分析 ======
moduleColors <- labels2colors(net$colors)
MEList = moduleEigengenes(dataExpr, colors = moduleColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

plotEigengeneNetworks(MEs,
                      "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)

table(moduleColors)

# 针对 Blue 模块的可视化
which.module="blue";
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(2,2.3,0.5,0.8))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

# TOM网络热图绘制
load(net$TOMFiles, verbose=T)
TOM <- as.matrix(TOM)
TOM[1:4,1:4]
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Call the plot function
table(moduleColors)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot, all genes",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))

# 输出 moduleLabels
Labels = as.data.frame(moduleLabels)
Colors = as.data.frame(moduleColors)
Labels = cbind(Labels,Colors)
Colors = cbind(id=rownames(Labels),Labels)
write.table(Colors,file="moduleLabels.txt",sep="\t",quote=F,row.names=F)
# ====================================================

# 特征基因关联临床特征 (保留B部分关于样本顺序一致性的优化)
cli = read.csv("Cli.csv", header = TRUE, check.names = FALSE)
rownames(cli) <- cli$ID
# 确保 MEs 与 cli 样本顺序一致
common_samples = intersect(rownames(MEs), rownames(cli))
MEs = MEs[common_samples, ]
cli = cli[common_samples, ]

moduleTraitCor <- cor(MEs, cli, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(common_samples))

pdf(file = "Module-trait.pdf", width = 7, height = 7)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
par(mar = c(10, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(cli),
               yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.65, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()
