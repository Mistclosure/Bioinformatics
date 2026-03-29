# ==============================================================================
# 1. 环境准备与数据读取
# ==============================================================================
setwd("/mnt/disk1/qiuzerui/downloads/Rloop/GSE146100")
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(devtools)
library(SingleR)
library(celldex) 
library(data.table)

rt_raw = fread("GSE146100_NormData.txt", data.table = FALSE)
rownames(rt_raw) = rt_raw[, 1]
rt = rt_raw[, -1]

# 创建 Seurat 对象
scRNA = CreateSeuratObject(rt, min.cells = 1, project = "W1", min.features = 1)

# 【核心修改】：由于 rt 是标准化后的，需手动将数据同步到 data 层
# Seurat v5 默认将输入存入 counts，而后续分析需要从 data 层提取标准化值
LayerData(scRNA, assay = "RNA", layer = "data") <- LayerData(scRNA, assay = "RNA", layer = "counts")

# ==============================================================================
# 2. 标准预处理流程 (替代 SCTransform)
# ==============================================================================
# 1. 寻找高变基因 (直接基于标准化的 data 层)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# 2. 数据缩放 (PCA 的必要前提)
# ScaleData 将每个基因的表达量进行中心化和标准化
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))

# ==============================================================================
# 3. 降维与聚类
# ==============================================================================
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

pc.num=1:30

scRNA <- RunTSNE(scRNA, dims=pc.num) %>% RunUMAP(dims=pc.num)

scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 0.5)

# ==============================================================================
# 4. 使用 SingleR 进行自动细胞类型鉴定
# ==============================================================================
refdata <- HumanPrimaryCellAtlasData()

# 【修改】：由于没有跑 SCT，此处 assay 应改为 "RNA"
testdata <- LayerData(scRNA, assay = "RNA", layer = "data")

clusters <- scRNA$seurat_clusters

cellpred <- SingleR(test = testdata, 
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels, 
                      stringsAsFactors = FALSE)

scRNA$singleR1 <- celltype[match(clusters, celltype$ClusterID), 'celltype']

# ==============================================================================
# 5. 可视化
# ==============================================================================
DimPlot(scRNA, group.by="singleR1", label=FALSE, label.size=5, reduction='tsne')
DimPlot(scRNA, group.by="orig.ident", label=FALSE, label.size=5, reduction='tsne')
