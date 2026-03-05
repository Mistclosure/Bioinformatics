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

load("scRNA.Rdata")
#计算Malignant cells
# 1. 加载所需的包
library(copykat)
library(dplyr)
library(stringr) # 用于拆分字符串

# 2. 提取全部细胞矩阵，并拼接细胞名与分组 ID
# 确保 Layers 已经合并，以提取完整的原始 counts 矩阵
scRNA <- JoinLayers(scRNA, assay = "RNA")
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

# 【核心修改】将细胞名和 id (如 P10, C3) 拼接成新的细胞名
# 使用 "___" (三个下划线) 作为分隔符，防止和原始细胞名中自带的单下划线冲突
# 假设此时 scRNA$orig.ident 已经映射为了 P1, P2 等格式
colnames(raw_counts) <- paste0(colnames(raw_counts), "___", scRNA$orig.ident)

# 3. 运行 CopyKAT
# 警告：现在是全量 7.4 万细胞运行，极度消耗内存和时间！建议务必在后台运行
copykat.res <- copykat(
  rawmat = as.matrix(raw_counts), 
  id.type = "S",           # "S" 代表基因符号 (Symbol)，即你矩阵里的基因名
  ngene.chr = 5,           # 每条染色体至少保留的基因数
  win.size = 25,           # 窗口大小
  KS.cut = 0.1,            # 判定非整倍体的严格程度，0.1 是保守推荐值
  sam.name = "lung_cancer", 
  distance = "euclidean", 
  n.cores = 16              # 请根据你服务器的 CPU 核心数进行调整，越大越快
)
save(copykat.res,file='copykat_res.Rdata')

# ==================== 修改开始 ====================
# 4. 提取预测结果并筛选恶性细胞 (非整倍体 Aneuploid)
pred <- as.data.frame(copykat.res$prediction)
malignant_joined_names <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 5. 利用正则表达式将拼接的细胞名分开成两列 (原始 barcode 和 分组 ID)
# 按照 "___" 分割，生成一个矩阵，第一列是原细胞名，第二列是 P1/C3 等 ID
split_names <- str_split(malignant_joined_names, "___", simplify = TRUE)
malig_barcodes <- split_names[, 1] # 真实的细胞 barcode
malig_ids <- split_names[, 2]      # 分组 ID (P1, M1...)

# --- 任务 A: 统计每个 ID 的恶性细胞数量，生成 Paint_Malignant cells.txt ---
malig_counts <- as.data.frame(table(malig_ids))
colnames(malig_counts) <- c("id", "number")

# 读取你的 Cli.csv (包含 ID 和 Original_Sample_ID)
cli_full <- read.csv("Cli.csv", header = TRUE)

# 通过 id (P1, M1) 将统计结果合并进表格，保证所有临床样本都在
meat <- merge(cli_full[, c("ID", "Original_Sample_ID")], malig_counts, 
              by.x = "ID", by.y = "id", all.x = TRUE)

# 把因没检测到恶性细胞而产生的 NA 替换为 0 (严谨补丁)
meat$number[is.na(meat$number)] <- 0

# 严格调整列顺序为：id, Original_Sample_ID, number (并将 ID 重命名为 id)
meat <- meat[, c("ID", "Original_Sample_ID", "number")]
colnames(meat)[1] <- "id"

# 保存 A 部分所需文件
write.table(meat, file = "Paint_Malignant cells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
print("恶性细胞统计文件 (Paint_Malignant cells.txt) 已生成完毕！")

# --- 任务 B: 提取全部恶性细胞，生成 C 部分所需的 Malignant cells.txt ---
# 根据 C 部分代码 `Malignant=read.table(..., row.names=1)` 和 `scRNA[,Malignant[,1]]`
# 说明该 txt 文件需要有行名，且第一数据列必须是细胞 barcode
malig_df <- data.frame(Barcode = malig_barcodes, row.names = malig_barcodes)

# 保存 C 部分所需文件
write.table(malig_df, file = "Malignant cells.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)
print("恶性细胞列表文件 (Malignant cells.txt) 已生成完毕！")
# ==================== 修改结束 ====================

#A----
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

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
scRNA$orig.ident <- factor(scRNA$orig.ident,levels=c("C1","C2","C3",
                                                     "M1","M2","M3",
                                                     "M4","M5","M6","M7","M8",
                                                     "M9","M10",
                                                     "P1","P2","P3","P4","P5",
                                                     "P6","P7","P8","P9","P10",
                                                     "P11","P12","P13","P14","P15",
                                                     "P16","P17"))
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
