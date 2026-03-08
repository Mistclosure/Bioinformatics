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
# 1. 加载所需的包
library(copykat)
library(dplyr)
library(stringr) # 用于拆分字符串

# 2. 提取全部细胞矩阵，并拼接细胞名与分组 ID
# 确保 Layers 已经合并，以提取完整的原始 counts 矩阵
scRNA <- JoinLayers(scRNA, assay = "RNA")
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

# 将细胞名和 id (如 P10, C3) 拼接成新的细胞名
colnames(raw_counts) <- paste0(colnames(raw_counts), "___", scRNA$orig.ident)


# ==================== 修改：按 orig.ident 分批运行与合并的核心代码 ====================

# 获取所有唯一的样本/批次 ID (如 C1, P1)
sample_ids <- unique(scRNA$orig.ident)

# 初始化一个空列表，用于存放各组的预测结果
pred_list <- list()

# 循环运行每一个批次/样本
for (i in seq_along(sample_ids)) {
  current_id <- sample_ids[i]
  print(paste0("========== 正在运行第 ", i, " 组 (共 ", length(sample_ids), " 组): ", current_id, " =========="))
  
  # 根据 orig.ident 提取当前批次的表达矩阵
  cell_idx <- which(scRNA$orig.ident == current_id)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]
  
  # 运行 CopyKAT
  sub_copykat_res <- copykat(
    rawmat = as.matrix(sub_counts), 
    id.type = "S",            
    ngene.chr = 5,            
    win.size = 25,            
    KS.cut = 0.1,             
    sam.name = paste0("copykat_", current_id), # 使用具体的 orig.ident 作为输出文件前缀
    distance = "euclidean", 
    n.cores = 32               
  )
  
  # 提取当前组的预测结果并存入列表
  # (加入了一个小的防错机制，以防个别细胞过少的样本跑不出结果而中断循环)
  if (!is.null(sub_copykat_res) && "prediction" %in% names(sub_copykat_res)) {
    pred_list[[current_id]] <- as.data.frame(sub_copykat_res$prediction)
  } else {
    print(paste0("警告：样本 ", current_id, " 未能返回有效预测结果，已跳过。"))
  }
  
  # 主动释放内存，防止后续循环爆内存
  rm(sub_counts, sub_copykat_res)
  gc()
}

print("========== 所有分组运行完毕，正在合并预测结果 ==========")

# 将所有组的预测结果按行合并为一个完整的数据框
pred <- do.call(rbind, pred_list)
# 保存合并后的预测结果以防万一
save(pred, file = 'copykat_merged_pred.Rdata')

# ==================== 修改结束 ====================


# ==================== 以下为未修改的原始代码 ====================
# 4. 提取预测结果并筛选恶性细胞 (非整倍体 Aneuploid)
# 注意：此时的 pred 已经是完整合并后的数据了，代码无缝衔接
malignant_joined_names <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 5. 利用正则表达式将拼接的细胞名分开成两列 (原始 barcode 和 分组 ID)
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
              by.x = "Original_Sample_ID", by.y = "id", all.x = TRUE)

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
malig_df <- data.frame(Barcode = malig_barcodes, row.names = malig_barcodes)

# 保存 C 部分所需文件
write.table(malig_df, file = "Malignant cells.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)
print("恶性细胞列表文件 (Malignant cells.txt) 已生成完毕！")

#A----
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# 1. 读取你整理好的 CSV 文件
cli_full = read.csv("Cli.csv", header = TRUE, check.names = FALSE)

##此步不要重复运行，否则会出错
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

meta$x <- factor(meta$orig.ident,levels=c("P1","P2","P3","P4","P5",
                                          "P6","P7","P8","P9","P10",
                                          "P11","P12","P13","P14","P15",
                                          "P16","P17","M1","M2","M3",
                                          "M4","M5","M6","M7","M8",
                                          "M9","M10","C1","C2","C3"))
ggplot(data = meta, aes(x = x, fill =singleR))+
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
