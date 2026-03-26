setwd("/mnt/disk1/qiuzerui/downloads/CRC/GSE178341/")
library(data.table)
library(Seurat)
library(stringr)
library(magrittr)
library(harmony) 
library(celldex)
data <- Read10X_h5("GSE178341_crc10x_full_c295v4_submit.h5")
scRNA = CreateSeuratObject(data, min.cells = 3, project =
                             "GSE178341",min.features =300)
# 计算线粒体比例 (这是优化的前提)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

### SCT 优化
# 1. method = "glmGamPoi" 大幅提速
# 2. vst.flavor = "v2" 提高准确性
# 3. vars.to.regress = "percent.mt" 消除线粒体干扰
scRNA <- SCTransform(scRNA, 
                     method = "glmGamPoi", 
                     vst.flavor = "v2", 
                     vars.to.regress = "percent.mt", 
                     verbose = FALSE)

### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

### Harmony
# 维持原样，SCT 配合 Harmony 效果很好
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", 
                    max.iter.harmony = 20, assay.use = "SCT")

pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num)

# FindNeighbors 增加 reduction="harmony" 显式指定，确保分群是基于矫正后的数据
scRNA <- FindNeighbors(scRNA, reduction="harmony", dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 1.5)
# 读取 CSV 文件
cluster <- read.csv("crc10x_full_c295v4_submit_cluster.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

scRNA@meta.data$cluster <- cluster$clMidwayPr[match(colnames(scRNA),cluster$sampleID)]
save(scRNA,'scRNA.Rdata')
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
    n.cores = 16               
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
cli_full <- read.csv("GSE178341_crc10x_full_c295v4_submit_metatables.csv", header = TRUE)

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

# ==================== 针对恶性细胞子集 (pbmc1) 的完整重跑流程 ====================

library(limma)
# 1. 读取恶性细胞列表并提取子集
Malignant = read.table("Malignant cells.txt", header=T, sep="\t", check.names=F, row.names=1)

# 确保从含有完整 metadata 的 scRNA 中提取
pbmc1 = scRNA[, rownames(Malignant)] 

# 2. 合并 Layers 并进行 SCTransform (重复全量数据的 SCT 逻辑)
pbmc1 <- JoinLayers(pbmc1)
# 建议：加入 method = "glmGamPoi" 提速，加入 vars.to.regress 回归线粒体影响
pbmc1 <- SCTransform(pbmc1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# 3. 准备 RNA 层 (专门用于你后续要跑的 CDSscore)
DefaultAssay(pbmc1) <- "RNA"
pbmc1 <- NormalizeData(pbmc1)
pbmc1 <- ScaleData(pbmc1)

# 4. 重新跑 PCA (基于子集的高变基因)
DefaultAssay(pbmc1) <- "SCT" # 降维必须回到 SCT
pbmc1 <- RunPCA(pbmc1, npcs = 50, verbose = FALSE)

# 5. 重新跑 Harmony (核心修改：纠正恶性细胞内部的病人差异)
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "orig.ident", 
                    assay.use = "SCT", max.iter.harmony = 20)

# 6. 重新跑降维与聚类 (基于 harmony 空间)
pc.num = 1:20 # 根据你之前设置的 20 个 dims
pbmc1 <- RunTSNE(pbmc1, reduction = "harmony", dims = pc.num) %>%
  RunUMAP(reduction = "harmony", dims = pc.num)

pbmc1 <- FindNeighbors(pbmc1, reduction = "harmony", dims = pc.num)
pbmc1 <- FindClusters(pbmc1, resolution = 1)

# 7. 可视化检查
# 此时建议分别按 cluster 和你之前的 MetastasisStatus 观察
p1 <- DimPlot(pbmc1, reduction = "umap", label = TRUE) + NoLegend()
p2 <- DimPlot(pbmc1, reduction = "umap", group.by = "orig.ident") # 检查 Harmony 效果

# 保存结果
save(pbmc1, file = 'Malignant.Rdata')

# ==================== 修改结束 ====================
