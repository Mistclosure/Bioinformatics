library(tidyverse)

# ================= 配置区域 (下次只改这里) =================
# 输入文件路径
input_file <- "counts/all_samples_featureCounts.txt" 
# 输出文件路径
output_file <- "Phf20_GSE82115_featurecounts.csv"
# ===========================================================

# 1. 读取数据
# comment.char = "#" 用于跳过 featureCounts 输出文件的第一行命令注释
# check.names = FALSE 极为重要！防止 R 自动把文件名里的特殊符号乱改
raw_counts <- read.table(input_file, header = TRUE, row.names = 1, 
                         comment.char = "#", check.names = FALSE)

# 2. 数据清洗与筛选
clean_counts <- raw_counts %>%
  # (A) 删除 featureCounts 自带的非 Count 列 (Chr, Start, End, Strand, Length)
  # 使用 select 排除法，这样比指定保留列更通用
  select(-any_of(c("Chr", "Start", "End", "Strand", "Length"))) %>%
  
  # (B) 清洗列名 (泛化核心)
  rename_with(.fn = function(x) {
    # 逻辑：STAR/featureCounts 的文件名通常包含固定的后缀和路径前缀
    # 我们用正则表达式把它们“替换为空”，剩下的就是样本名
    
    x %>%
      # 1. 去掉 .Aligned 及其后面所有的内容 (比如 .Aligned.sortedByCoord.out.bam)
      str_remove("\\.Aligned.*") %>%
      
      # 2. 去掉 alignments. 及其前面所有的内容 (包括 X.home... 等路径乱码)
      # 这里假设你的文件名里有 "alignments" 这个词。如果没有，可以用 "/" 
      str_remove(".*alignments/") %>%
      
      # 3. (可选) 如果还有残留的 X. 前缀 (R 读取数字开头文件时自动加的)，去掉它
      str_remove("^X\\.")
  })

# 3. 检查结果 (在控制台打印前几行和列名看看对不对)
cat("清洗后的列名预览:\n")
print(colnames(clean_counts))

# 4. 保存文件
write.csv(clean_counts, output_file, quote = FALSE)

cat("\n处理完成！文件已保存为:", output_file)
