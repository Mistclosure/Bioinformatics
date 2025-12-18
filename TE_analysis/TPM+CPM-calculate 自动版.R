# ==============================================================================
# 综合处理脚本：双输出模式 (单独输出 CPM 文件 和 TPM 文件)
# ==============================================================================

library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# ==============================
# 1. 参数与路径设置
# ==============================
counts_file <- "Phf20_GSE_counts.csv"

# 注释文件路径
gene_gtf_path <- "\\\\wsl.localhost\\Ubuntu\\home\\qiuzerui\\annotationMv38\\gencode.vM38.annotation_PRI.gtf"
te_gtf_path   <- "\\\\wsl.localhost\\Ubuntu\\home\\qiuzerui\\annotationMv38\\m39_TE.gtf"

# 定义两个输出文件名
output_cpm_file <- "Phf20_GSE82115_CPM.csv"  # 输出1: 包含 Counts 和 CPM
output_tpm_file <- "Phf20_GSE82115_TPM.csv"  # 输出2: 包含 Length 和 TPM
samplename <- c('shNT_rep1','shNT_rep2','shNT_rep3','shPHF20_rep1','shPHF20_rep2','shPHF20_rep3')

# ==============================
# 2. 读取数据与准备注释 (Symbol转换)
# ==============================
message(paste0("[", Sys.time(), "] 正在读取 Counts 文件..."))
counts_df <- fread(counts_file)
colnames(counts_df)[2:ncol(counts_df)] <- samplename
# --- 2.1 加载 GTF 提取 Gene Symbol (为了让 CPM 文件也有 Symbol) ---
message(paste0("[", Sys.time(), "] 正在加载 GTF 以匹配 Gene Symbol..."))
gene_gtf <- import(gene_gtf_path)
gene_map <- unique(as.data.frame(mcols(gene_gtf)[, c("gene_id", "gene_name")]))

# --- 2.2 ID 转换: ENSG -> Symbol ---
# 只有 Gene (ENSG开头) 需要转，TE (带冒号) 保持原样
counts_df <- left_join(counts_df, gene_map, by = c("RepeatID" = "gene_id"))

# 如果匹配到了 gene_name，就用 gene_name 替换 RepeatID
counts_df$RepeatID <- ifelse(
  !is.na(counts_df$gene_name), 
  counts_df$gene_name, 
  counts_df$RepeatID
)
counts_df$gene_name <- NULL # 删除临时列

message("   - ID 转换完成 (ENSG -> Gene Symbol)")

# ==============================
# 3. 计算并输出 CPM (独立文件)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算 CPM..."))

# 复制一份数据用于计算 CPM
df_cpm <- copy(counts_df)

# 识别样本列 (排除 ID 列)
sample_cols <- setdiff(names(df_cpm), "RepeatID")

# 计算 CPM
for (col in sample_cols) {
  # 库大小 = 该样本所有 Gene + TE 的 reads 总和
  library_size <- sum(df_cpm[[col]], na.rm = TRUE)
  
  # 计算 CPM
  cpm_val <- (df_cpm[[col]] / library_size) * 1e6
  
  # 添加新列 (保留2位小数)
  df_cpm[[paste0(col, "_CPM")]] <- round(cpm_val, 2)
}

# --- 输出文件 1: CPM ---
message(paste(">>> 正在导出 CPM 文件:", output_cpm_file))

# 可选：调整列顺序，让 Count 和 CPM 挨在一起
# cols_order_cpm <- c("RepeatID")
# for (sample in sample_cols) {
#   cols_order_cpm <- c(cols_order_cpm, sample, paste0(sample, "_CPM"))
# }
# df_cpm <- df_cpm %>% select(any_of(cols_order_cpm))

write.csv(df_cpm, output_cpm_file, row.names = FALSE)
message("   ✅ CPM 文件已保存 (包含 Raw Counts 和 CPM)")

# ==============================
# 4. 计算长度 (Gene + TE) - 为 TPM 做准备
# ==============================
message(paste0("[", Sys.time(), "] 正在计算基因长度 (用于 TPM)..."))

# --- 4.1 Gene 长度 ---
gene_exons <- gene_gtf[gene_gtf$type == "exon"]
gene_exons_list <- split(gene_exons, mcols(gene_exons)$gene_id)
gene_widths <- sum(width(reduce(gene_exons_list)))
gene_len_df <- data.frame(id = names(gene_widths), length = as.numeric(gene_widths))
# 再次匹配 Symbol (因为 ID 已经是 Symbol 了，我们需要用 gene_id 关联一下或者直接用 map)
# 为了准确，我们这里还是用 gene_id 关联 map 再转成 symbol
gene_len_df <- left_join(gene_len_df, gene_map, by = c("id" = "gene_id"))
gene_len_df$final_id <- ifelse(!is.na(gene_len_df$gene_name), gene_len_df$gene_name, gene_len_df$id)
gene_len_final <- gene_len_df[, c("final_id", "length")]

# --- 4.2 TE 长度 ---
te_gtf <- import(te_gtf_path)
te_mcols <- mcols(te_gtf)
# 构造与 Count 表一致的 TE ID
te_ids <- paste(te_mcols$gene_id, te_mcols$family_id, te_mcols$class_id, sep = ":")
mcols(te_gtf)$te_unique_id <- te_ids
te_list <- split(te_gtf, mcols(te_gtf)$te_unique_id)
te_widths <- sum(width(reduce(te_list)))
te_len_final <- data.frame(final_id = names(te_widths), length = as.numeric(te_widths))

# --- 4.3 合并长度表 ---
all_lengths <- rbind(gene_len_final, te_len_final)

# ==============================
# 5. 计算并输出 TPM (独立文件)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算 TPM..."))

# 将长度信息合并回主表
# 注意：counts_df 里的 RepeatID 已经是 Symbol/TE_ID 了，可以直接与 all_lengths$final_id 匹配
df_tpm <- left_join(counts_df, all_lengths, by = c("RepeatID" = "final_id"))

# 移除无长度的行
df_tpm <- df_tpm[!is.na(df_tpm$length), ]

# 计算 TPM
calculate_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  scaling_factor <- sum(rpk) / 1e6
  return(rpk / scaling_factor)
}

# 仅对样本列计算 (排除 ID 和 length)
tpm_calc_cols <- setdiff(names(df_tpm), c("RepeatID", "length"))

# 创建一个新的只包含 TPM 的表 (如果你想保留Counts也可以，这里演示保留 Counts+TPM)
# 如果只想保留 TPM，可以新建一个 dataframe
for (col in tpm_calc_cols) {
  tpm_val <- calculate_tpm(df_tpm[[col]], df_tpm$length)
  df_tpm[[paste0(col, "_TPM")]] <- round(tpm_val, 4)
}
#调整length列顺序
cols_order_tpm <- c("RepeatID",'length')
tpm_calc_cols <- setdiff(names(df_tpm), c("RepeatID", "length"))
for (sample in tpm_calc_cols) {
  cols_order_tpm <- c(cols_order_tpm, sample)
}
df_tpm <- df_tpm %>% select(any_of(cols_order_tpm))
#--- 输出文件 2: TPM ---
message(paste(">>> 正在导出 TPM 文件:", output_tpm_file))

# (可选)这里我们只保留 ID, Length 和 TPM 列 (去除 Raw Counts 以保持文件纯净，按需调整)
# 也就是：RepeatID, length, Sample1_TPM, Sample2_TPM...
# cols_keep_tpm <- c("RepeatID", "length", grep("_TPM$", names(df_tpm), value = TRUE))
# df_tpm_out <- df_tpm %>% select(all_of(cols_keep_tpm))

write.csv(df_tpm, output_tpm_file, row.names = FALSE)
message("   ✅ TPM 文件已保存 (包含 Length 和 TPM)")

message("========================================================")
message("🎉 全部完成！已生成两个独立文件：")
message(paste("1.", output_cpm_file))
message(paste("2.", output_tpm_file))
message("========================================================")
message("========================================================")
