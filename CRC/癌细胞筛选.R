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
# 读取 CSV 文件
cluster <- read.csv("crc10x_full_c295v4_submit_cluster.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

scRNA@meta.data$cluster <- cluster$clMidwayPr[match(colnames(scRNA),cluster$sampleID)]
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


# ==================== 修改：按 orig.ident 分为两大部分运行与合并的核心代码 ====================

# 获取所有唯一的样本/批次 ID
sample_ids <- unique(scRNA$orig.ident)

# 将批次尽量平均分为两大部分
mid_idx <- ceiling(length(sample_ids) / 2)
batch1_ids <- sample_ids[1:mid_idx]
batch2_ids <- sample_ids[(mid_idx + 1):length(sample_ids)]

# 将两部分存入列表以便循环
split_batches <- list(Part1 = batch1_ids, Part2 = batch2_ids)

# 初始化一个空列表，用于存放两大部分的预测结果
pred_list <- list()

# 循环运行这两大批次
for (i in seq_along(split_batches)) {
  current_batch_name <- names(split_batches)[i]
  current_batch_ids <- split_batches[[i]]
  
  print(paste0("========== 正在运行 ", current_batch_name, " (共2部分) =========="))
  print(paste0("当前部分包含的样本: ", paste(current_batch_ids, collapse=", ")))
  
  # 提取当前部分的细胞 barcode 索引
  cell_idx <- which(scRNA$orig.ident %in% current_batch_ids)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]
  
  # 运行 CopyKAT
  sub_copykat_res <- copykat(
    rawmat = as.matrix(sub_counts), 
    id.type = "S",            
    ngene.chr = 5,            
    win.size = 25,            
    KS.cut = 0.1,             
    sam.name = paste0("copykat_", current_batch_name), # 输出前缀改为 Part1/Part2
    distance = "euclidean", 
    n.cores = 16               
  )
  
  # 提取当前大部分的预测结果并存入列表
  if (!is.null(sub_copykat_res) && "prediction" %in% names(sub_copykat_res)) {
    pred_list[[current_batch_name]] <- as.data.frame(sub_copykat_res$prediction)
  } else {
    print(paste0("警告：", current_batch_name, " 未能返回有效预测结果，已跳过。"))
  }
  
  # 主动释放内存，清理当前大部分的稠密矩阵和运行结果
  rm(sub_counts, sub_copykat_res)
  gc()
}

print("========== 两大部分运行完毕，正在合并预测结果 ==========")

# 将两部分的预测结果按行合并为一个完整的数据框
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
