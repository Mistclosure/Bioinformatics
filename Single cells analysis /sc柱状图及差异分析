# ==============================================================================
# 0. 准备工作：提取 PBMC 对象
# ==============================================================================

pbmc <- PBMC

# 确保 Metadata 中的 Group 是因子，并设置合理的顺序 (方便画图)
# 顺序：25度 -> 30度 -> 4度 (或者按你喜欢的逻辑)
pbmc$Group <- factor(pbmc$Group, levels = c("RT_25C", "TN_30C", "Cold_4C"))

# ==============================================================================
# 1. 绘制细胞类型百分比堆积柱状图 (带百分比标签)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(scales) # 用于显示百分比标签

print("📊 正在绘制细胞比例堆积柱状图 (带标签)...")

# 1.1 构建绘图数据
# 统计每个 Group 中每种 cell_type 的细胞数量
cell_stats <- pbmc@meta.data %>%
  group_by(Group, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(percentage = count / sum(count)) # 计算小数比例

# --- 【核心修改】准备标签文本 ---
# 使用 scales::percent 将小数转换为百分比字符串，accuracy = 0.1 保留一位小数
cell_stats$label_text <- scales::percent(cell_stats$percentage, accuracy = 0.1)

# 1.2 绘图
p_bar <- ggplot(cell_stats, aes(x = Group, y = percentage, fill = cell_type)) +
  # 绘制柱状图层
  geom_bar(stat = "identity", position = "fill", width = 0.7) + 
  
  # --- 【核心修改】添加文字标签图层 ---
  # data = subset(...) 用于过滤掉太小的色块(例如小于3%)，避免文字重叠挤在一起
  # 修改为 > 0.01 (即显示大于 1% 的标签)
  geom_text(data = subset(cell_stats, percentage > 0.01), 
            aes(label = label_text),
            position = position_fill(vjust = 0.5),
            size = 2.5, # 建议同时把字体改小一点，防止挤不下
            color = "black") +
  
  # 设置Y轴和标签
  scale_y_continuous(labels = scales::percent) + 
  labs(title = "PBMC Cell Type Proportion by Group",
       x = "Condition", 
       y = "Percentage of Cells",
       fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

# 1.3 保存
ggsave(filename = "PBMC_CellType_Proportion_Labeled.png", plot = p_bar, width = 10, height = 7, path = data_dir)
print(paste("   ✅ 带标签的比例图已保存: PBMC_CellType_Proportion_Labeled.png"))

# ==============================================================================
# 2. Monocyte 差异表达分析 (Differential Expression)
# ==============================================================================
print("🔬 正在进行 Monocyte 差异分析...")

# 2.1 锁定 Monocytes 细胞
# ScType 注释的结果里可能包含 "Classical Monocytes", "Non-classical Monocytes" 等
# 我们先看看有哪些包含 "Monocyte" 的类型
monocyte_types <- grep("Monocyte", unique(pbmc$cell_type), value = TRUE)

if (length(monocyte_types) == 0) {
  stop("⚠️ 未在 PBMC 中找到名称包含 'Monocyte' 的细胞类型，请检查 cell_type 注释结果！")
}

print(paste("   已选定分析的细胞类型:", paste(monocyte_types, collapse = ", ")))

# 提取这些细胞的子集
pbmc_mono <- subset(pbmc, subset = cell_type %in% monocyte_types)

# 切换 Identity 为 Group，方便做组间比较
Idents(pbmc_mono) <- "Group"

# 再次确保图层合并 (防止 FindMarkers 报错)
pbmc_mono <- JoinLayers(pbmc_mono)


# --- 定义一个通用的差异分析与绘图函数 ---
run_deg_analysis <- function(obj, ident.1, ident.2, label_prefix) {
  
  print(paste("   >>> 正在对比:", ident.1, "vs", ident.2))
  
  # 1. 运行 FindMarkers
  # logfc.threshold = 0.25 (默认值，筛选差异明显的)
  # min.pct = 0.1 (只保留至少在10%细胞中表达的基因)
  deg_table <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2, 
                           test.use = "wilcox", verbose = FALSE)
  
  # 2. 添加基因名列 (方便保存)
  deg_table$gene <- rownames(deg_table)
  
  # 3. 筛选显著差异基因 (p_val_adj < 0.05)
  sig_genes <- deg_table %>% filter(p_val_adj < 0.05)
  
  # 区分上调和下调 (基于 ident.1)
  # avg_log2FC > 0 表示在 ident.1 (Cold_4C) 中高表达
  up_genes <- sig_genes %>% filter(avg_log2FC > 0)
  down_genes <- sig_genes %>% filter(avg_log2FC < 0)
  
  print(paste("       显著上调基因数:", nrow(up_genes)))
  print(paste("       显著下调基因数:", nrow(down_genes)))
  print(paste("       Top 5 上调基因:", paste(head(up_genes$gene, 5), collapse=",")))
  
  # 4. 保存结果到 CSV
  file_name <- paste0(data_dir, "/DEG_Monocytes_", label_prefix, ".csv")
  write.csv(deg_table, file = file_name, row.names = FALSE)
  print(paste("       结果已保存至:", file_name))
  
  # 5. 绘制简易火山图 (Volcano Plot)
  # 标记显著性分类
  deg_table$diff <- "NO"
  deg_table$diff[deg_table$avg_log2FC > 0.5 & deg_table$p_val_adj < 0.05] <- "UP"
  deg_table$diff[deg_table$avg_log2FC < -0.5 & deg_table$p_val_adj < 0.05] <- "DOWN"
  
  # 提取 Top 10 显著基因用于标注
  top_genes <- rbind(
    head(subset(deg_table, diff == "UP"), 10),
    head(subset(deg_table, diff == "DOWN"), 10)
  )
  
  library(ggrepel)
  
  p_vol <- ggplot(deg_table, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = diff), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_text_repel(data = top_genes, aes(label = gene), max.overlaps = 20) +
    labs(title = paste("Volcano: Monocytes -", label_prefix),
         subtitle = paste(ident.1, "vs", ident.2),
         x = "log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0("Volcano_Monocytes_", label_prefix, ".png"), 
         plot = p_vol, width = 7, height = 6, path = data_dir)
}

# --- 执行两次比较 ---

# 场景 1: 4度(寒冷) vs 25度(常温)
# 注意：ident.1 是实验组(分子)，ident.2 是对照组(分母)
# 结果中的 LogFC > 0 代表在 Cold_4C 中上调
tryCatch({
  run_deg_analysis(pbmc_mono, ident.1 = "Cold_4C", ident.2 = "RT_25C", label_prefix = "4C_vs_25C")
}, error = function(e) { message("❌ 4C vs 25C 分析出错: ", e$message) })

# 场景 2: 4度(寒冷) vs 30度(热中性)
tryCatch({
  run_deg_analysis(pbmc_mono, ident.1 = "Cold_4C", ident.2 = "TN_30C", label_prefix = "4C_vs_30C")
}, error = function(e) { message("❌ 4C vs 30C 分析出错: ", e$message) })

print("🎉 分析全部完成！请查看生成的 .csv 表格和 .png 图片。")
