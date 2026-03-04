#fig1----
#data----
setwd("/mnt/disk1/qiuzerui/downloads/Rloop/GSE131907")
library(data.table)
library(Seurat)
library(stringr)
library(magrittr)
library(harmony) 
library(celldex)
UMI = fread("GSE131907_Lung_Cancer_raw_UMI_matrix.txt",data.table = F)
rownames(UMI) = UMI[,1]
UMI = UMI[,-1]
scRNA = CreateSeuratObject(UMI, min.cells = 3, project =
                             "GSE131907",min.features =300)
ann = read.table("GSE131907_Lung_Cancer_cell_annotation.txt", header=T,
                 sep="\t", check.names=F, row.names=1)
ann = ann[colnames(scRNA),]
scRNA@meta.data$Sample = ann$Sample

# ==================== 新增修改：统一 ID ====================
# 将 GSE131907 的真实样本名覆盖掉默认的 project name，使其与 GSE123904 格式对齐
scRNA@meta.data$orig.ident = scRNA@meta.data$Sample
# ===========================================================

meta = scRNA@meta.data
id = meta[meta$Sample=="LUNG_T09" | meta$Sample=="LUNG_T08" |
            meta$Sample=="LUNG_T25" | meta$Sample=="LUNG_T06" |
            meta$Sample=="LUNG_T34" | meta$Sample=="LUNG_T31" |
            meta$Sample=="LUNG_T19" | meta$Sample=="LUNG_T20" |
            meta$Sample=="LUNG_T18" | meta$Sample=="LUNG_T28" |
            meta$Sample=="NS_06" | meta$Sample=="NS_16" |
            meta$Sample=="NS_02" | meta$Sample=="NS_19" |
            meta$Sample=="NS_17" | meta$Sample=="NS_04" |
            meta$Sample=="NS_13" | meta$Sample=="NS_03" |
            meta$Sample=="NS_12" | meta$Sample=="NS_07" |
            meta$Sample=="LUNG_T30",]
scRNA1 = scRNA[,rownames(id)]
save(scRNA1,file='scRNA1.Rdata')


setwd("/mnt/disk1/qiuzerui/downloads/Rloop/GSE123904/GSE123904_RAW/")
samples=list.files("./")
samples
dir <- file.path('./',samples)
names(dir) <- samples
names(dir)=str_split(names(dir),'_',simplify = T)[,3]
scRNAlist <- list()
for(i in 1:length(dir)){
  A = fread(dir[i])
  A = as.data.frame(A)
  rownames(A) = A[,1]
  A = as.data.frame(t(A[,-1]))
  scRNAlist[[i]] <- CreateSeuratObject(A, min.cells = 3, project =
                                         names(dir)[i],min.features =300)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]],
                                    scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],
                                    scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]],scRNAlist[[9]]))
scRNA = merge(scRNA1,y=scRNA2)
### SCT
scRNA <- SCTransform(scRNA)
### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony =
                      20,assay.use = "SCT")
pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 1.5)
table(scRNA@meta.data$seurat_clusters)
library(SingleR)
library(celldex)

# 对于肺癌，通常推荐使用 Human Primary Cell Atlas (HPCA) 或 BlueprintEncode
# 1. 下载参考数据集
refdata <- HumanPrimaryCellAtlasData()
# 1. 定义 testdata：在 v5 中使用 'layer' 取代 'slot'
# 注意：对于 SCTransform 处理后的数据，归一化矩阵存储在 'data' 图层
testdata <- GetAssayData(scRNA, assay = "SCT", layer = "data")

# 2. 定义 clusters：获取聚类结果 (这一步语法与 v4 保持一致)
# 但建议显式指定，确保获取的是最新的聚类列
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels =
                      refdata$label.main,
                     clusters = clusters,
                    assay.type.test = "logcounts", assay.type.ref =
                      "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred),
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
setwd('/mnt/disk1/qiuzerui/downloads/Rloop')
save(scRNA,file='scRNA.Rdata')

#计算Malignant cells
# 1. 加载所需的包
library(copykat)
library(dplyr)

# 2. 提取上皮细胞亚集 (显著降低计算量并提高准确率)
# 注意：HPCA 参考库里上皮细胞的标签通常是 "Epithelial_cells"
# 你可以先用 table(scRNA$singleR) 确认一下具体的拼写
epi_cells <- subset(scRNA, subset = singleR == "Epithelial_cells")

# 确保 Layers 已经合并，以提取完整的原始 counts 矩阵
epi_cells <- JoinLayers(epi_cells, assay = "RNA")
raw_counts <- GetAssayData(epi_cells, assay = "RNA", layer = "counts")

# 3. 运行 CopyKAT
# 警告：这一步依然需要一段时间（可能几个小时），建议在 tmux 或 screen 后台运行
copykat.res <- copykat(
  rawmat = as.matrix(raw_counts), 
  id.type = "S",           # "S" 代表基因符号 (Symbol)，即你矩阵里的基因名
  ngene.chr = 5,           # 每条染色体至少保留的基因数
  win.size = 25,           # 窗口大小
  KS.cut = 0.1,            # 判定非整倍体的严格程度，0.1 是保守推荐值
  sam.name = "lung_cancer", 
  distance = "euclidean", 
  n.cores = 64              # 请根据你服务器的 CPU 核心数进行调整，越大越快
)

# 4. 提取预测结果并筛选恶性细胞 (非整倍体 Aneuploid)
pred <- as.data.frame(copykat.res$prediction)
# 获取被鉴定为恶性细胞的 barcode (细胞名)
malignant_barcodes <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 5. 统计每个样本的恶性细胞数量
# 从 Seurat 对象中提取这些恶性细胞的元数据
malignant_meta <- scRNA@meta.data[malignant_barcodes, ]

# 统计每个新 ID (P1, P2, M1...) 的细胞数量
malig_counts <- as.data.frame(table(malignant_meta$orig.ident))
colnames(malig_counts) <- c("id", "number")

# 6. 拼装符合原作者代码要求的表格格式
# 原代码要求：第一列是简化 ID(P1), 第二列是原始 ID(LUNG_T09), 第三列是数目(number)
# 我们读取你之前整理的 Cli.csv 来获取原始 ID
cli_full <- read.csv("Cli.csv", header = TRUE)

# 将原始 ID 合并进统计结果表
meat <- merge(malig_counts, cli_full[, c("ID", "Original_Sample_ID")], 
              by.x = "id", by.y = "ID", all.x = TRUE)

# 严格调整列顺序为：id, Original_Sample_ID, number
meat <- meat[, c("id", "Original_Sample_ID", "number")]

# 7. 保存为 txt 文件
write.table(meat, file = "Paint_Malignant cells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

print("恶性细胞统计文件已生成完毕！")


#A----
library(ComplexHeatmap)
library(circlize)

# 1. 读取你整理好的 CSV 文件
cli_full = read.csv("Cli.csv", header = TRUE, check.names = FALSE)

# 2. 进行原始 ID 到新 ID 的映射
# 构建一个以 Original_Sample_ID 为名，ID (P1, M1...) 为值的字典向量
id_map <- setNames(cli_full$ID, cli_full$Original_Sample_ID)
# 映射替换 scRNA 对象的 orig.ident
scRNA$orig.ident <- unname(id_map[scRNA$orig.ident])

# 3. 整理供后续热图使用的 cli 格式
cli <- cli_full
rownames(cli) <- cli$ID
# 删掉不需要画在热图上的前三列：ID, Original_Sample_ID, Dummy
cli <- cli[, -c(1, 2, 3)] 
# （重要补救）你的 CSV 表头是 "Origins"，但原代码要求 "Sample Origins"，这里用代码帮你强行对齐，以保证后续代码不报错
colnames(cli)[colnames(cli) == "Origins"] <- "Sample Origins" 

value = rnorm(19)
colnames(cli)

# ==================== 以下为原作者代码，一字未改 ====================

ha = HeatmapAnnotation(df = cli,
                       col = list(
                         `Data Source` = c("GSE131907" = "#E69394" ,
                                           "GSE123904" = "#BEBADA"),
                         `Sample Origins` = c("Primary" = "#B3E2CD" ,
                                              "Distant Metastasis" = "#E4D4B7", "Chemotherapy" = "#ECCFC0"),
                         Smoking = c("Never smoker" = "#2DB600","Current\\nsmoker" = "#EDB48E", "Former smoker" = "#E6E600"),
                         EGFR = c("WT" = "#7FC97F", "Mut" = "#FDC086",
                                  "na" = "#A9B7B7"),
                         KRAS = c("WT" = "#7FC97F", "Mut" = "#FDC086",
                                  "na" = "#A9B7B7"),
                         TP53 = c("WT" = "#7FC97F", "Mut" = "#FDC086",
                                  "na" = "#A9B7B7"),
                         Stage = c("Stage I" = "#EFF3FF", "Stage II" =
                                     "#B8D4E6", "Stage III" = "#64A9D3", "Stage IV" = "#2A7AB7")
                       ))
draw(ha)

meat=read.table("Paint_Malignant cells.txt", header=T, sep="\t",
                check.names=F)
meat$x <- factor(meat$id,levels=c("P1","P2","P3","P4","P5",
                                  "P6","P7","P8","P9","P10",
                                  "P11","P12","P13","P14","P15",
                                  "P16","P17","M1","M2","M3",
                                  "M4","M5","M6","M7","M8",
                                  "M9","M10","C1","C2","C3"))
meat$number = log2(meat$number+1)

ggplot(meat, aes(x=x, y=number, group = 1)) +
  geom_line(size=2,color="#CBD5E8")+
  geom_point(size=3)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(x="",y="log2(Malignant Cell Number)")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

meta = scRNA@meta.data
meta = meta[,c(1,10)]
for (i in 1:nrow(meat)) {
  x = meat[i,1]
  y = meat[i,2]
  meta$orig.ident[which(meta$orig.ident == y)] <- x
}

meta$x <- factor(meta$orig.ident,levels=c("P1","P2","P3","P4","P5",
                                          "P6","P7","P8","P9","P10",
                                          "P11","P12","P13","P14","P15",
                                          "P16","P17","M1","M2","M3",
                                          "M4","M5","M6","M7","M8",
                                          "M9","M10","C1","C2","C3"))
ggplot(data = meta, aes(x = x, fill =singleR1))+
  geom_bar(stat = 'count',position = 'fill')+labs(y = "Cell\nProportion(%)" , x="")+
  scale_fill_manual(values = c( "#80B1D3","#BC80BD" , "#FB8072"
                                ,"#8DD3C7", "#FFFFB3",
                                "#FDB462" ,"#D9D9D9","#FCCDE5",
                                "#BABADA","#B3DE69"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

#B----可运行
DimPlot(scRNA, group.by="singleR", label.size=5, reduction='umap')
DimPlot(scRNA, group.by="orig.ident", label.size=5, reduction='umap')


#C----
library(limma)
Malignant=read.table("Malignant cells.txt", header=T, sep="\t",
                     check.names=F, row.names=1)
pbmc1 = scRNA[,Malignant[,1]]
Count = as.data.frame(pbmc1@assays[["RNA"]]@counts)
meta = pbmc1@meta.data
pbmc1 <- CreateSeuratObject(counts = Count)
pbmc1@meta.data$Type = B$orig.ident
### SCT
pbmc1 <- SCTransform(pbmc1)
#PCA
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1))
pbmc1 <- RunTSNE(pbmc1, dims=pc.num) %>% RunUMAP(dims=1:20)
pbmc1 <- FindNeighbors(pbmc1, dims = 1:20)
pbmc1 <- FindClusters(pbmc1, resolution = 1)
DimPlot(pbmc1, reduction = "umap")
save(pbmc1,file='Malignant.Rdata')
