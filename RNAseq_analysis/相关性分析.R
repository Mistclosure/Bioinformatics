setwd('E:/Documents-3/TE通路富集/')
# ==============================================================================
# 1. 环境准备与包安装
# ==============================================================================
# 检查常用 CRAN 包
cran_packages <- c("tidyverse", "data.table", "ggrepel", "readxl")
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran)

# 检查并安装 Bioconductor 包 (GSVA)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 确保安装的是较新版本的 GSVA
if (!require("GSVA", quietly = TRUE)) {
  cat("正在安装 GSVA 包...\n")
  BiocManager::install("GSVA")
}

library(tidyverse)
library(data.table)
library(ggrepel)
library(readxl)
library(GSVA)  # 加载 ssGSEA 核心包

# 检查 GSVA 版本，如果是旧版本提示升级
if (packageVersion("GSVA") < "1.50.0") {
  warning("您的 GSVA 版本较旧 (< 1.50.0)。建议使用 BiocManager::install('GSVA') 升级以支持新语法。")
}

# ==============================================================================
# 2. 读取数据
# ==============================================================================
cat("Step 1: 正在读取数据文件...\n")

# A. 读取 hERV 分组表
herv_meta_file <- "JCI121476.sdt3.xlsx" 
if (file.exists(herv_meta_file)) {
  herv_meta <- read_excel(herv_meta_file, sheet = 1, .name_repair = "minimal")
} else {
  stop("找不到 JCI121476.sdt3.xlsx，请检查文件名。")
}

# B. 读取 hERV 表达矩阵 (Table 12 UQN)
herv_exp_file <- "JCI121476.sdt12.txt"
if(file.exists(herv_exp_file)){
  herv_exp <- fread(herv_exp_file)
} else {
  stop("找不到 JCI121476.sdt12.txt")
}

# C. 读取基因表达矩阵 (HiSeqV2)
gene_exp_file <- "HiSeqV2" 
if(file.exists(gene_exp_file)){
  gene_exp_raw <- fread(gene_exp_file, header = TRUE)
} else {
  stop("找不到 HiSeqV2 文件")
}

# ==============================================================================
# 3. 数据清洗与匹配
# ==============================================================================
cat("Step 2: 正在清洗数据与匹配样本 ID...\n")

# --- hERV 矩阵处理 ---
if (!"Sample_ID" %in% colnames(herv_exp)) colnames(herv_exp)[1] <- "Sample_ID"
herv_df <- herv_exp %>% column_to_rownames("Sample_ID")

# --- 基因矩阵处理 (转置) ---
gene_names <- gene_exp_raw[[1]] 
gene_exp_mat <- t(as.matrix(gene_exp_raw[, -1]))
colnames(gene_exp_mat) <- gene_names
gene_df <- as.data.frame(gene_exp_mat)

# --- ID 统一 ---
clean_ids <- function(ids) {
  ids <- gsub("[_.]", "-", ids)
  substr(ids, 1, 15)
}
rownames(herv_df) <- clean_ids(rownames(herv_df))
rownames(gene_df) <- clean_ids(rownames(gene_df))

# --- 取交集 ---
common_samples <- intersect(rownames(herv_df), rownames(gene_df))
cat(paste0("匹配成功：共找到 ", length(common_samples), " 个共有样本。\n"))

herv_clean <- herv_df[common_samples, ]
gene_clean <- gene_df[common_samples, ]

# --- 基因数据 Log 检查 ---
max_val <- max(gene_clean, na.rm = TRUE)
if (max_val > 50) {
  cat("检测到线性数据，执行 log2(x+1) 转换...\n")
  gene_clean <- log2(gene_clean + 1)
} else {
  cat("检测到 Log2 数据，无需转换。\n")
}

# ==============================================================================
# 4. 计算 ssGSEA Scores (全新 GSVA 语法)
# ==============================================================================
cat("Step 3: 使用新版 GSVA 计算 ssGSEA Scores...\n")

col_g1 <- "RIG-I-like up (Group 1)"
col_g2 <- "RIG-I-like down (Group 2)"

# 1. 提取 ID 列表
g1_ids <- herv_meta %>% filter(!!sym(col_g1) == TRUE | !!sym(col_g1) == "True") %>% pull(1)
g2_ids <- herv_meta %>% filter(!!sym(col_g2) == TRUE | !!sym(col_g2) == "True") %>% pull(1)

# 仅保留矩阵中存在的 ID
g1_ids <- intersect(g1_ids, colnames(herv_clean))
g2_ids <- intersect(g2_ids, colnames(herv_clean))

cat(paste0("Group 1 hERV 数量: ", length(g1_ids), "\n"))
cat(paste0("Group 2 hERV 数量: ", length(g2_ids), "\n"))

# 2. 构建 Gene Sets 列表
herv_gene_sets <- list(
  Group1_ssGSEA = g1_ids,
  Group2_ssGSEA = g2_ids
)

# 3. 准备矩阵 (Features x Samples)
# GSVA 要求输入矩阵为：行=基因/特征，列=样本
herv_matrix_t <- t(as.matrix(herv_clean))

# 4. 运行 ssGSEA (New Syntax)
cat("正在运行 gsva()... \n")

if (packageVersion("GSVA") >= "1.50.0") {
  # --- [新版语法] ---
  # 先构建参数对象 ssgseaParam
  # alpha: 权重指数，默认 0.25 是 ssGSEA 的标准配置
  # normalize: 是否标准化结果，默认 TRUE
  ssgsea_params <- ssgseaParam(exprData = herv_matrix_t, 
                               geneSets = herv_gene_sets,
                               alpha = 0.25,
                               normalize = TRUE)
  
  # 运行计算
  ssgsea_res <- gsva(ssgsea_params)
  
} else {
  # --- [兼容旧版] ---
  cat("检测到旧版 GSVA，使用旧语法...\n")
  ssgsea_res <- gsva(herv_matrix_t, 
                     herv_gene_sets, 
                     method = "ssgsea", 
                     kcdf = "Gaussian")
}

# 转置结果 (变成 行=样本，列=Score)
scores <- as.data.frame(t(ssgsea_res))
# 此时 scores 包含 Group1_ssGSEA 和 Group2_ssGSEA 两列

# ==============================================================================
# 5. 全基因组相关性分析
# ==============================================================================
cat("Step 4: 计算全基因组相关性 (Based on ssGSEA Scores)...\n")

target_gene <- "PHF20"
if (!target_gene %in% colnames(gene_clean)) {
  match <- grep("PHF20", colnames(gene_clean), value = TRUE)
  if(length(match) > 0) target_gene <- match[1]
}

# 定义计算函数
calc_cor_safe <- function(gene_df, score_vec, label) {
  res <- apply(gene_df, 2, function(x) {
    test <- cor.test(as.numeric(x), score_vec, method = "pearson")
    c(Correlation = unname(test$estimate), P_value = test$p.value)
  })
  res_df <- as.data.frame(t(res))
  # 调整列顺序，Gene 放首位
  res_df <- res_df %>% mutate(Gene = rownames(.)) %>% select(Gene, Correlation, P_value)
  res_df$Target <- label
  return(res_df)
}

# --- 计算 Group 1 ---
cat("正在计算 Group 1 相关性...\n")
res_g1 <- calc_cor_safe(gene_clean, scores$Group1_ssGSEA, "Group 1 (ssGSEA)")

# --- 计算 Group 2 ---
cat("正在计算 Group 2 相关性...\n")
res_g2 <- calc_cor_safe(gene_clean, scores$Group2_ssGSEA, "Group 2 (ssGSEA)")

# 保存结果
write.csv(res_g1, "hERV_Group1_ssGSEA_Correlation.csv", row.names = FALSE)
write.csv(res_g2, "hERV_Group2_ssGSEA_Correlation.csv", row.names = FALSE)
cat("表格已保存至 hERV_Group1/2_ssGSEA_Correlation.csv\n")

# ==============================================================================
# 6. PHF20 排名与绘图
# ==============================================================================
cat("\nStep 5: PHF20 分析与绘图...\n")

# 获取排名
get_rank <- function(df) {
  df %>% arrange(desc(Correlation)) %>% mutate(Rank = 1:n()) %>% filter(Gene == target_gene)
}

r1 <- get_rank(res_g1)
r2 <- get_rank(res_g2)

cat(paste0("PHF20 在 Group 1 中的排名: ", r1$Rank, " / ", nrow(res_g1), " (R = ", round(r1$Correlation,3), ")\n"))
cat(paste0("PHF20 在 Group 2 中的排名: ", r2$Rank, " / ", nrow(res_g2), " (R = ", round(r2$Correlation,3), ")\n"))

# 绘图函数
plot_rank <- function(res_df, group_name, output_filename) {
  plot_data <- res_df %>% 
    arrange(desc(Correlation)) %>% 
    mutate(Rank = 1:n()) %>% 
    mutate(Highlight = ifelse(Gene == target_gene, "Target", "Other"))
  
  target_point <- plot_data %>% filter(Highlight == "Target")
  
  p <- ggplot(plot_data, aes(x = Rank, y = Correlation, color = Highlight)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_point(data = target_point, color = "red", size = 4) +
    geom_text_repel(data = target_point, 
                    label = paste0(target_gene, "\nR=", round(target_point$Correlation, 3), "\nRank=", target_point$Rank),
                    box.padding = 1.5, nudge_y = 0.2, nudge_x = 1000, color = "black") +
    scale_color_manual(values = c("Other" = "grey70", "Target" = "red")) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = paste0("Gene Correlation Rank vs ", group_name),
         subtitle = "Based on New GSVA (ssGSEA)",
         x = "Rank (Positive -> Negative)", y = "Pearson Correlation")
  
  ggsave(output_filename, plot = p, width = 8, height = 6, dpi = 300)
}

plot_rank(res_g1, "Group 1 (ssGSEA)", "Plot_Group1_ssGSEA_Rank.png")
plot_rank(res_g2, "Group 2 (ssGSEA)", "Plot_Group2_ssGSEA_Rank.png")

cat("图片已保存: Plot_Group1_ssGSEA_Rank.png, Plot_Group2_ssGSEA_Rank.png\n")
