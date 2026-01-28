# ==============================================================================
# TE 差异分析脚本 (最终修复版)
# 功能：
# 1. 筛选 TE
# 2. 拆分 SubFamily, Family, Class
# 3. 计算 Log2FC 和 PValue
# 4. 仅保留分类、统计值和 CPM 数据 (去除原始 RepeatID, Symbol, Counts)
# ==============================================================================

# 加载必要的包
library(data.table)
library(dplyr)
library(tidyr)
library(stringr) # 确保加载 stringr

# --- 设置输入输出文件 ---
input_cpm <- "Phf20_1.23_CPM.csv"
output_te_stats <- "Phf20_1.23_TE_Stats_Clean.csv"

message(paste0("[", Sys.time(), "] 正在读取 CPM 文件..."))
df <- fread(input_cpm)

# 1. 筛选 TE (至少包含2个冒号)
#    确保使用 stringr 包的函数
df_te <- df[stringr::str_count(df$RepeatID, ":") >= 2, ]
message(paste0("   - 筛选出 TE 数量: ", nrow(df_te)))

# 2. 拆分 ID 为三列 (SubFamily, Family, Class)
#    注意：separate 后转为 data.frame 避免 data.table 的索引兼容性问题
df_te <- df_te %>%
  separate(RepeatID, c("SubFamily", "Family", "Class"), sep = ":", remove = FALSE, extra = "merge") %>%
  as.data.frame()

# 3. 定义 CPM 列名 (用于计算和输出)
ctrl_cpm_cols <- c("L1MKL2609676-Scr_1_Mixt_CPM", "L1MKL2609677-Scr_2_Mixt_CPM")
treat_cpm_cols <- c("L1MKL2609678-Phf20_1_Mixt_CPM", "L1MKL2609679-Phf20_2_Mixt_CPM")
all_cpm_cols <- c(ctrl_cpm_cols, treat_cpm_cols)

# 检查列名是否存在
if (!all(all_cpm_cols %in% names(df_te))) {
  stop("错误：CPM 列名匹配失败，请检查文件中的列名！")
}

# ==============================
# 4. 计算 Log2FC (向量化计算)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算统计指标..."))

# 提取 CPM 数据矩阵
mat_ctrl <- df_te[, ctrl_cpm_cols]
mat_treat <- df_te[, treat_cpm_cols]

# 计算均值 (添加 0.01 防止 log(0))
pseudo <- 0.01
mean_ctrl <- rowMeans(mat_ctrl, na.rm = TRUE)
mean_treat <- rowMeans(mat_treat, na.rm = TRUE)

# 计算 Log2FC
df_te$Log2FC <- round(log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo), 4)

# ==============================================================================
# 5. 计算 P-value (包含函数定义与执行)
# ==============================================================================
message(paste0("[", Sys.time(), "] 正在执行 T-test 计算..."))

# --- A. 定义计算函数 (工人) ---
# 该函数接收一行数据，提取指定列，进行 Student's t-test (var.equal=TRUE)
calc_pval <- function(x, c_cols, t_cols) {
  # 强制转换为数值，防止数据框格式导致的计算错误
  v_c <- as.numeric(x[c_cols])
  v_t <- as.numeric(x[t_cols])
  
  # 安全检查：如果两组数据完全一样（方差均为0），t.test 会报错
  # 这种情况直接返回 P=1 (无显著差异)
  if (var(v_c, na.rm = TRUE) == 0 && var(v_t, na.rm = TRUE) == 0) return(1)
  
  # 执行检验
  tryCatch({
    # var.equal = TRUE 指定使用普通 Student's t-test
    res <- t.test(v_t, v_c, var.equal = TRUE)
    return(res$p.value)
  }, error = function(e) {
    # 如果发生其他不可预见的错误，返回 NA
    return(NA_real_)
  })
}

# --- B. 初始化与准备 ---
# 先创建空列，确保后续 select() 不会因为找不到列而崩溃
df_te$PValue <- NA_real_

# 确保输入 apply 的是纯粹的数据框
df_te_mat <- as.data.frame(df_te)

# --- C. 调用 apply 执行计算 (流水线) ---
# 1 代表按行处理
df_te$PValue <- apply(df_te_mat, 1, calc_pval, 
                      c_cols = ctrl_cpm_cols, 
                      t_cols = treat_cpm_cols)

# --- D. 后期处理 ---
# 转换为数值并保留 5 位小数
df_te$PValue <- round(as.numeric(df_te$PValue), 5)

message("   ✅ P-value 计算完成！")
# ==============================
# 6. 整理输出 (移除 RepeatID, Symbol, Counts)
# ==============================

# 定义最终要保留的列：分类 + 统计 + CPM
final_cols <- c(
  "SubFamily", "Family", "Class", 
  "Log2FC", "PValue", 
  all_cpm_cols
)

# 筛选列
df_final <- df_te %>% select(all_of(final_cols))

# 导出文件
message(paste(">>> 正在导出精简版统计文件:", output_te_stats))
write.csv(df_final, output_te_stats, row.names = FALSE)

message("   ✅ 全部完成！已生成文件：", output_te_stats)
