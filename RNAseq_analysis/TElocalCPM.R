# ==============================================================================
# Phf20 TElocal 专用 CPM 计算脚本
# ==============================================================================
setwd('/mnt/windowsdata/qiuzerui/Phf20-26.1.23/')
library(data.table)
library(dplyr)
library(rtracklayer)

# ==============================
# 1. 路径设置
# ==============================
input_file  <- "Phf20_1.23_TElocal_locus_counts.csv"
output_file <- "Phf20_1.23_TElocal_CPM_Calculated.csv"

# 基因注释路径（用于将 ENSG 转换成 Symbol）
gene_gtf_path <- "/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38/gencode.vM38.annotation_PRI.gtf"

# ==============================
# 2. 读取数据与注释映射
# ==============================
message(paste0("[", Sys.time(), "] 正在读取 Counts 文件..."))
counts_df <- fread(input_file)

message(paste0("[", Sys.time(), "] 正在加载 GTF 以匹配 Gene Symbol..."))
gene_gtf <- import(gene_gtf_path)
# 提取 ID 和 Symbol 的对应关系
gene_map <- unique(as.data.frame(mcols(gene_gtf)[, c("gene_id", "gene_name")]))

# 关联 Symbol：如果是基因则显示名字，如果是 TE 则显示原始 ID
counts_df <- left_join(counts_df, gene_map, by = c("RepeatID" = "gene_id"))
counts_df$Symbol <- ifelse(!is.na(counts_df$gene_name), counts_df$gene_name, counts_df$RepeatID)

# 整理列顺序：RepeatID, Symbol, 然后是原始 Counts
counts_df$gene_name <- NULL
col_order <- c("RepeatID", "Symbol", setdiff(names(counts_df), c("RepeatID", "Symbol")))
df_final <- counts_df %>% select(all_of(col_order))

# ==============================
# 3. 核心计算：CPM (修复数值类型报错)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算 CPM..."))

# 自动识别样本列
sample_cols <- setdiff(names(df_final), c("RepeatID", "Symbol"))

# 循环每一列进行 CPM 计算
for (col in sample_cols) {
  message("   - 正在处理样本: ", col)
  
  # --- 关键修复：强制转换列为数值型 ---
  # 使用 as.numeric 转换，防止 'character' 参数无效报错
  counts_vec <- as.numeric(df_final[[col]])
  
  # 如果转换过程中产生了 NA (说明原数据有非数字内容)，将其填为 0
  counts_vec[is.na(counts_vec)] <- 0
  
  # 获取该样本的总 Read 数
  total_counts <- sum(counts_vec, na.rm = TRUE)
  
  # 如果总和为 0，防止除以零报错
  if (total_counts == 0) {
      message("     ⚠️ 警告: 样本 ", col, " 总 Count 为 0，跳过 CPM 计算。")
      df_final[[paste0(col, "_CPM")]] <- 0
      next
  }
  
  # 计算 CPM 公式: (count / total_count) * 1,000,000
  cpm_values <- (counts_vec / total_counts) * 1e6
  
  # 将结果存入新列，列名为：原名_CPM
  df_final[[paste0(col, "_CPM")]] <- round(cpm_values, 2)
}
# ==============================
# 4. 保存结果
# ==============================
message(paste0("[", Sys.time(), "] >>> 正在导出结果至: ", output_file))
write.csv(df_final, output_file, row.names = FALSE)

message("========================================================")
message("✅ CPM 计算完成！")
message("结果包含：")
message("1. 原始 Counts 列")
message("2. 新增的 _CPM 后缀列")
message("3. 基因 Symbol 注释列")
message("========================================================")
