# ==============================================================================
# 第一步：环境准备与安装包
# ==============================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 整合两组代码需要的包
required_packages <- c("GSVA", "limma", "ggplot2", "dplyr", "data.table", "org.Hs.eg.db", "AnnotationDbi")
for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE)){
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 第二步：加载 DEG 数据并构建 Signature
# ==============================================================================
# 请确保路径正确
# setwd('E:/Documents-3/力学signature/') 

# 1. 读取您的差异表达分析结果文件
deg_file <- "顺铂敏感和抗性上皮细胞株DEGs_results.csv"
res <- read.csv(deg_file, stringsAsFactors = FALSE)

# 2. 定义力学相关关键词
mech_keywords <- c("ACT", "MYO", "TUB", "KRT", "VIM", "DES", "NEFL", "MAP", 
                   "CDH", "ITG", "CLDN", "OCLN", "JAM", "CTN", 
                   "COL", "LAMA", "FN1", "FBN", "FBLN", "MMP", "TIMP", "LOX", "TNC", 
                   "TAGLN", "CNN", "FLNA", "VCL", "TLN", "ZYX", "PXN")
pattern <- paste(mech_keywords, collapse = "|")

# 3. 筛选 Signature 基因
signature_candidates <- res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2 & baseMean > 100) %>%
  filter(grepl(pattern, gene_symbol, ignore.case = TRUE)) %>%
  arrange(desc(log2FoldChange))

# 4. 提取最终 Signature (Top 15)
target_genes <- head(signature_candidates$gene_symbol, 15)

print("--- 最终使用的 Mechanics Signature ---")
print(target_genes)

# 构建 Gene Set 列表
gene_sets <- list(Mechanics_Score = target_genes)

# ==============================================================================
# 第三步 (修复版)：读取 TCGA-BLCA 数据并进行稳健的 ID 转换
# ==============================================================================

library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

print("正在读取 TCGA 表达矩阵...")
# 1. 读取数据
blca_file <- "TCGA-BLCA.star_tpm.tsv"
tpm_data <- fread(blca_file, data.table = FALSE)

# 2. 提取 Ensembl ID 并去除版本号
# 假设第一列是 ID (通常列名是 Ensembl_ID 或 gene_id)
# 您的报错显示第一列名可能是 "Ensembl_ID"，或者是无名列
raw_ids <- tpm_data[, 1]
clean_ids <- sub("\\..*", "", raw_ids)

# 3. ID 映射 (Ensembl -> Symbol)
print("正在将 Ensembl ID 转换为 Gene Symbol...")
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = clean_ids, 
                                  columns = "SYMBOL", 
                                  keytype = "ENSEMBL")

# 4. 构建合并用的数据框
# 移除原始的第一列 ID，只保留数值数据(样本列)
expr_values <- tpm_data[, -1] 
# 将清洗后的 ID 加回去用于匹配
expr_values$ensembl_clean <- clean_ids

# 5. 合并 Symbol
# inner join (只保留能匹配到 Symbol 的行)
merged_df <- merge(expr_values, gene_map, by.x = "ensembl_clean", by.y = "ENSEMBL", all = FALSE)

# 移除没有 Symbol 的行 (NA)
merged_df <- merged_df[!is.na(merged_df$SYMBOL), ]

# 6. 智能聚合 (核心修复步骤)
print("正在合并重复的 Gene Symbol (这步很关键)...")

# 找出所有的数值列 (即样本列)
# 排除 ensembl_clean, SYMBOL, ENSEMBL 等非数值列
numeric_cols <- sapply(merged_df, is.numeric)
# 确保只对数值列进行聚合，且以 SYMBOL 为分组依据
# 使用 aggregate 的公式接口，但为了稳健，我们先构建一个清洁的子集
data_for_agg <- merged_df[, numeric_cols]
data_for_agg$SYMBOL <- merged_df$SYMBOL

# 执行聚合 (加入 na.rm = TRUE 防止任何潜在的 NA 扩散)
expr_clean <- aggregate(. ~ SYMBOL, data = data_for_agg, FUN = mean, na.rm = TRUE)

# 7. 重建最终矩阵
# 第一列现在是 SYMBOL，把它变成行名
expr_mat <- as.matrix(expr_clean[, -1])
rownames(expr_mat) <- expr_clean$SYMBOL

# 8. 最终清洗：处理矩阵中可能残留的 NA
# 如果整行全是 NA，删除；如果有零星 NA，补 0 (防止 GSVA 报错)
if(any(is.na(expr_mat))) {
  print("警告：矩阵中仍检测到 NA 值，正在进行清洗...")
  # 替换 NA 为 0 (通常 log 后的 0 代表不表达)
  expr_mat[is.na(expr_mat)] <- 0
}

print("转换后的矩阵预览 (确保没有 Ensembl_ID 列，且数值正常)：")
print(expr_mat[1:5, 1:5])

# ==============================================================================
# 第四步：Log2 转化与运行 ssGSEA
# ==============================================================================

# Log2 转化 (使用 ID 转换后的 expr_mat)
if(max(expr_mat, na.rm=TRUE) > 50) {
  print("检测到原始 TPM 值，正在进行 log2(TPM + 1) 转化...")
  expr_mat <- log2(expr_mat + 1)
} else {
  print("数据似乎已经 log 转化过，直接使用。")
}

library(GSVA)

# 简单验证 Signature 匹配情况
overlap <- intersect(rownames(expr_mat), gene_sets$Mechanics_Score)
print(paste("成功匹配到的 Signature 基因数量:", length(overlap)))

if(length(overlap) < 5) {
  stop("错误：匹配到的基因太少！请检查 ID 转换是否成功。")
}

print("正在构建 ssGSEA 参数对象...")

# 使用新版 ssgseaParam (去除 kcdf)
params <- ssgseaParam(exprData = expr_mat, 
                      geneSets = gene_sets, 
                      minSize = 1,        
                      maxSize = Inf,
                      alpha = 0.25,      
                      normalize = TRUE)

print("正在计算分数...")
ssgsea_score <- gsva(params)

# 结果整理
scores_df <- as.data.frame(t(ssgsea_score))
scores_df$SampleID <- rownames(scores_df)

print("计算完成！")
print(head(scores_df))

# ==============================================================================
# 第五步：预测与验证 (Correlation Analysis)
# ==============================================================================

# 提取分数向量
mech_vec <- scores_df$Mechanics_Score

# 验证的基因列表
target_check_genes <- c("CD70", "CD274", "IL33", "TGFBI", "TIGIT", "CTLA4")

# 确保基因在矩阵中
valid_genes <- intersect(target_check_genes, rownames(expr_mat))

if(length(valid_genes) > 0){
  print("--- 目标基因相关性验证结果 ---")
  
  for(g in valid_genes){
    gene_exp <- expr_mat[g, ]
    test <- cor.test(mech_vec, gene_exp, method = "pearson")
    
    cat(sprintf("Gene: %s | Correlation: %.3f | P-value: %.3e\n", 
                g, test$estimate, test$p.value))
    
    if(g == "CD70"){
      plot_df <- data.frame(Score = mech_vec, Expression = gene_exp)
      p <- ggplot(plot_df, aes(x = Score, y = Expression)) +
        geom_point(alpha = 0.5, color = "darkblue") +
        geom_smooth(method = "lm", color = "red") +
        theme_minimal() +
        labs(title = "Mechanics Score vs CD70 Expression",
             x = "Mechanics Score (ssGSEA)", y = "CD70 Expression (log2 TPM)")
      print(p)
      ggsave("Mechanics_vs_CD70_Correlation.png", width = 6, height = 5)
    }
  }
} else {
  print("警告：在表达矩阵中未找到目标验证基因。")
}

# 全基因组扫描
print("正在进行全基因组关联扫描 (Top 5000 变异基因)...")
top_var_genes <- names(sort(apply(expr_mat, 1, var), decreasing = TRUE))[1:5000]
cor_res <- data.frame(Gene = character(), Cor = numeric(), Pval = numeric(), stringsAsFactors = FALSE)

pb <- txtProgressBar(min = 0, max = length(top_var_genes), style = 3)
for(i in 1:length(top_var_genes)){
  g <- top_var_genes[i]
  if(g %in% target_genes) next 
  
  test <- cor.test(mech_vec, expr_mat[g, ], method = "pearson")
  cor_res[i, ] <- list(g, test$estimate, test$p.value)
  setTxtProgressBar(pb, i)
}
close(pb)

# 展示 Top 正相关基因
top_predicted <- cor_res %>% 
  filter(Pval < 0.05) %>% 
  arrange(desc(Cor)) %>% 
  head(20)

print("--- 预测的 Top 20 未知关联基因 ---")
print(top_predicted)
write.csv(top_predicted, "Predicted_Mechanics_Related_Genes.csv", row.names = FALSE)
