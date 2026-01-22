# ==============================================================================
# 1. 加载必要的包
# ==============================================================================
packages <- c("tidyverse", "data.table", "ggrepel", "readxl")
if (any(!packages %in% rownames(installed.packages()))) {
  install.packages(packages[!packages %in% rownames(installed.packages())])
}

library(tidyverse)
library(data.table)
library(ggrepel)
library(readxl)

# ==============================================================================
# 2. 读取数据
# ==============================================================================
cat("Step 1: 正在读取数据文件...\n")

# A. 读取 hERV 分组表
herv_meta_file <- "JCI121476.sdt3.xlsx" # 请确认文件名
if (file.exists(herv_meta_file)) {
  # 读取 Excel，保留列名原始格式
  herv_meta <- read_excel(herv_meta_file, sheet = 1, .name_repair = "minimal")
} else {
  stop("找不到 JCI121476.sdt3.xlsx，请检查文件名。")
}

# B. 读取 hERV 表达矩阵 (Table 12)
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

# --- ID 统一与清洗 ---
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

# --- [关键] 检查基因数据是否需要 Log2 ---
max_val <- max(gene_clean, na.rm = TRUE)
if (max_val > 50) {
  cat("检测到线性数据 (Max > 50)，正在执行 log2(x+1) 转换...\n")
  gene_clean <- log2(gene_clean + 1)
} else {
  cat("检测到 Log2 数据 (Max < 50)，无需转换。\n")
}

# ==============================================================================
# 4. 计算 Group 1 和 Group 2 评分
# ==============================================================================
cat("Step 3: 计算 hERV Group Scores...\n")

col_g1 <- "RIG-I-like up (Group 1)"
col_g2 <- "RIG-I-like down (Group 2)"

# 提取 ID
g1_ids <- herv_meta %>% filter(!!sym(col_g1) == TRUE | !!sym(col_g1) == "True") %>% pull(1)
g2_ids <- herv_meta %>% filter(!!sym(col_g2) == TRUE | !!sym(col_g2) == "True") %>% pull(1)

# 仅保留存在的 ID
g1_ids <- intersect(g1_ids, colnames(herv_clean))
g2_ids <- intersect(g2_ids, colnames(herv_clean))

# 计算 Score
scores <- data.frame(
  G1_Score = rowMeans(herv_clean[, g1_ids, drop=FALSE], na.rm = TRUE),
  G2_Score = rowMeans(herv_clean[, g2_ids, drop=FALSE], na.rm = TRUE)
)

# ==============================================================================
# 5. 全基因组相关性分析 (分开输出)
# ==============================================================================
cat("Step 4: 计算全基因组相关性...\n")

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
  res_df$Gene <- rownames(res_df)
  res_df$Target <- label
  return(res_df)
}

# --- 计算 Group 1 ---
res_g1 <- calc_cor_safe(gene_clean, scores$G1_Score, "Group 1 (Immunosuppressive)")
# [修改点] 将 Gene 列放到第一位
res_g1 <- res_g1 %>% select(Gene, Correlation, P_value, Target)

# --- 计算 Group 2 ---
res_g2 <- calc_cor_safe(gene_clean, scores$G2_Score, "Group 2 (Immunostimulatory)")
# [修改点] 将 Gene 列放到第一位
res_g2 <- res_g2 %>% select(Gene, Correlation, P_value, Target)

# --- [修改点] 分别保存表格 ---
write.csv(res_g1, "hERV_Group1_Correlation.csv", row.names = FALSE)
write.csv(res_g2, "hERV_Group2_Correlation.csv", row.names = FALSE)

cat("结果已保存:\n1. hERV_Group1_Correlation.csv\n2. hERV_Group2_Correlation.csv\n")

# ==============================================================================
# 6. PHF20 排名与绘图
# ==============================================================================
cat("\nStep 5: PHF20 分析与绘图...\n")

# 打印 PHF20 信息
p1 <- res_g1 %>% filter(Gene == target_gene)
p2 <- res_g2 %>% filter(Gene == target_gene)

# 计算排名
r1 <- res_g1 %>% arrange(desc(Correlation)) %>% mutate(Rank=1:n()) %>% filter(Gene == target_gene)
r2 <- res_g2 %>% arrange(desc(Correlation)) %>% mutate(Rank=1:n()) %>% filter(Gene == target_gene)

cat(paste0("PHF20 Group 1 排名: ", r1$Rank, " / ", nrow(res_g1), " (R=", round(r1$Correlation,3), ")\n"))
cat(paste0("PHF20 Group 2 排名: ", r2$Rank, " / ", nrow(res_g2), " (R=", round(r2$Correlation,3), ")\n"))

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
         x = "Rank (Positive -> Negative)", y = "Pearson Correlation")
  
  ggsave(output_filename, plot = p, width = 8, height = 6, dpi = 300)
}

plot_rank(res_g1, "Group 1", "Plot_Group1_Rank.png")
plot_rank(res_g2, "Group 2", "Plot_Group2_Rank.png")

cat("图片已保存: Plot_Group1_Rank.png, Plot_Group2_Rank.png\n")
