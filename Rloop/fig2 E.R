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
library(celldex) # SingleR 引用数据库所在的包

# 读取标准化后的表达矩阵
rt = read.table("GSE146100_NormData.txt", header=T, sep="\t", check.names=F, row.names=1)

# 创建 Seurat 对象 (Seurat v5 兼容 v4 的创建方式)
scRNA = CreateSeuratObject(rt, min.cells = 3, project = "W1", min.features = 300)

# ==============================================================================
# 2. 标准预处理流程 (SCTransform)
# ==============================================================================
# 使用 SCTransform 进行标准化、寻找高变基因并回归干扰因素
# 在 Seurat v5 中，SCTransform 依然是处理单细胞数据的首选
scRNA <- SCTransform(scRNA, verbose = FALSE)

# ==============================================================================
# 3. 降维与聚类
# ==============================================================================
# 运行 PCA 降维
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

# 设置用于后续分析的维度数
pc.num=1:30

# 运行 tSNE 和 UMAP 可视化算法
scRNA <- RunTSNE(scRNA, dims=pc.num) %>% RunUMAP(dims=pc.num)

# 构建 SNN 图并寻找细胞群 (此处修正了原代码的连字错误)
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 0.5)

# ==============================================================================
# 4. 使用 SingleR 进行自动细胞类型鉴定
# ==============================================================================
# 获取参考数据集 (人主细胞图谱)
refdata <- HumanPrimaryCellAtlasData()

# 【Seurat v5 修改】: 使用 LayerData 提取数据层
# 原 v4 语法: GetAssayData(scRNA, slot="data") 
# v5 推荐使用 layer="data" 来访问经过标准化后的矩阵
testdata <- LayerData(scRNA, assay = "SCT", layer = "data")

# 提取当前的聚类信息
clusters <- scRNA$seurat_clusters

# 运行 SingleR (基于 Cluster 模式进行鉴定)
cellpred <- SingleR(test = testdata, 
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    # 注意：SingleR 内部默认处理 log 后的数据
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

# 整理鉴定结果
celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels, 
                      stringsAsFactors = FALSE)

# 将鉴定结果映射回 Seurat 对象的 metadata 中
# 在 v5 中，直接使用 scRNA$name 赋值比操作 @meta.data 更简洁
scRNA$singleR1 <- celltype[match(clusters, celltype$ClusterID), 'celltype']

# ==============================================================================
# 5. 可视化
# ==============================================================================
# 按 SingleR 鉴定结果绘制 tSNE 图
DimPlot(scRNA, group.by="singleR1", label=FALSE, label.size=5, reduction='tsne')

# 按样本来源绘制 tSNE 图
DimPlot(scRNA, group.by="orig.ident", label=FALSE, label.size=5, reduction='tsne')
