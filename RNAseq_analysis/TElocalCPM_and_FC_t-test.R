# ==============================================================================
# Phf20 TElocal 专用 CPM 计算脚本
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/downloads/Phf8_GSE212779')
library(data.table)
library(dplyr)
library(rtracklayer)

# ==============================
# 1. 路径设置
# ==============================
input_file  <- "Phf8_GSE212779_TElocal_locus_counts.csv"
output_file <- "Phf8_GSE212779_TElocal_locus_CPM.csv"

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
# ==============================================================================
# TElocal 差异分析脚本 (灵活样本量与顺序版)
# 格式要求：RepeatID 为 locus:subfamily:family:class
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# --- 1. 设置路径与手动分组变量 ---
input_cpm       <- "Phf8_GSE212779_TElocal_locus_CPM.csv"
output_te_stats <- "Phf8_GSE212779_TElocal_locus_CPM_FC.csv"

# 【核心配置】：在此处直接列出所有样本的完整列名
# 脚本会自动识别样本量（不论是 2x2, 3x3 还是更多）
ctrl_samples  <- c("SRR21106101_CPM", "SRR21106102_CPM", "SRR21106103_CPM")
treat_samples <- c("SRR21106098_CPM", "SRR21106099_CPM", "SRR21106100_CPM")

# --- 2. 读取数据与验证 ---
message(paste0("[", Sys.time(), "] 正在读取数据并验证样本名..."))
df <- fread(input_cpm)

# 合并所有需要计算的列
all_specified_cols <- c(ctrl_samples, treat_samples)

# 自动检查列名是否存在，防止拼写错误
missing_cols <- setdiff(all_specified_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste("❌ 错误：在 CSV 中找不到以下列名，请检查拼写或后缀：\n", 
             paste(missing_cols, collapse = ", ")))
}

# --- 3. 筛选并拆分 TE ID (4层架构) ---
# TElocal locus 模式包含 3 个冒号 (locus:sub:fam:class)
df_te <- df[stringr::str_count(df$RepeatID, ":") == 3, ]
message(paste0("   - 筛选出有效 TE 数量: ", nrow(df_te)))

df_te <- df_te %>%
  separate(RepeatID, 
           into = c("Locus", "SubFamily", "Family", "Class"), 
           sep = ":", 
           remove = FALSE, 
           extra = "merge") %>%
  as.data.frame()

# --- 4. 数据类型安全转换 ---
# 强制将选定列转为数值型，防止计算报错
for(col in all_specified_cols) {
  df_te[[col]] <- as.numeric(df_te[[col]])
}
df_te[all_specified_cols][is.na(df_te[all_specified_cols])] <- 0

# --- 5. 计算 Log2FC ---
# 公式：$$Log_2FC = \log_2(\frac{\text{Mean}_{Treat} + \epsilon}{\text{Mean}_{Ctrl} + \epsilon})$$
message(paste0("[", Sys.time(), "] 正在计算 Log2FC..."))
pseudo <- 0.01
mean_ctrl  <- rowMeans(as.matrix(df_te[, ctrl_samples]), na.rm = TRUE)
mean_treat <- rowMeans(as.matrix(df_te[, treat_samples]), na.rm = TRUE)

df_te$Log2FC <- round(log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo), 4)

# --- 6. 计算 P-value (T-test) ---
message(paste0("[", Sys.time(), "] 正在执行 T-test 计算..."))

calc_pval_te <- function(x, c_cols, t_cols) {
  v_c <- as.numeric(x[c_cols])
  v_t <- as.numeric(x[t_cols])
  
  # 如果两组内部方差均为0（数值完全相同），直接返回 1 避免 T-test 报错
  if (var(v_c) == 0 && var(v_t) == 0) return(1)
  
  tryCatch({
    # 使用标准 Student's t-test
    res <- t.test(v_t, v_c, var.equal = TRUE)
    return(res$p.value)
  }, error = function(e) return(NA_real_))
}

# 动态传入识别到的列名
df_te$PValue <- apply(df_te, 1, calc_pval_te, 
                      c_cols = ctrl_samples, 
                      t_cols = treat_samples)

df_te$PValue <- round(as.numeric(df_te$PValue), 5)

# --- 7. 整理输出结果 ---
# 保留四层分类、统计指标及原始 CPM 数据
final_cols <- c("Locus", "SubFamily", "Family", "Class", "Log2FC", "PValue", all_specified_cols)
df_final <- df_te %>% select(all_of(final_cols))

write.csv(df_final, output_te_stats, row.names = FALSE)
message("✅ 全部完成！已处理 ", length(ctrl_samples), " 个对照 vs ", length(treat_samples), " 个实验样本。")
