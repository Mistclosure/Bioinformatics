# ==============================================================================
# 适配 Kat8 P60 WT vs KO 数据的完整流程
# 特性：自动处理生物学重复整合 (Harmony) + 智能文件名解析
# ==============================================================================
set.seed(42)
# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx) 
library(Matrix)    
library(scales)
library(harmony) # 强烈建议安装：install.packages("harmony") 用于整合生物学重复
library(HGNChelper) 
library(dplyr)
# --- 加载去双胞专用包 ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")

library(scDblFinder)
library(SingleCellExperiment)

# 加载 ScType 核心函数
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# ------------------------------------------------------------------------------
# 1. 设置路径与样本信息 (根据你的截图修改)
# ------------------------------------------------------------------------------
# !!! 修改为你存放这4个文件夹的路径 !!!
data_dir <- "D:\\qiuzerui\\kat8p60-2" 
setwd(data_dir)
# 根据截图中的文件夹名定义
sample_folders <- c("kat8-P60-WT-1", "kat8-P60-WT-2", "kat8-P60-Y90C-KO-1", "kat8-P60-Y90C-KO-2")

# ------------------------------------------------------------------------------
# 2. 读取数据、解析文件名元数据
# ------------------------------------------------------------------------------
sc_list <- list()

print("🚀 步骤1/6: 开始读取数据并解析实验分组...")

for (folder in sample_folders) {
  full_path <- file.path(data_dir, folder)
  
  if(!dir.exists(full_path)) {
    message(paste("⚠️ 跳过：未找到文件夹", folder)); next
  }
  
  print(paste("   正在处理:", folder))
  
  tryCatch({
    # 读取 10X 数据
    counts <- Read10X(data.dir = full_path)
    
    # 处理可能的 list 结构
    if (is.list(counts) && !is(counts, "dgCMatrix")) {
      counts <- if ("Gene Expression" %in% names(counts)) counts$`Gene Expression` else counts[[1]]
    }
    
    # 创建 Seurat 对象
    # 注意：project 参数直接用文件夹名，方便后续追踪
    sc_obj <- CreateSeuratObject(counts = counts, project = folder, min.cells = 3, min.features = 200)
    
    # --- 核心：智能解析文件名 ---
    # 格式示例: kat8-P60-WT-1
    # 格式示例: kat8-P60-Y90C-KO-1
    
    # 判断分组 (Group)
    if (grepl("WT", folder)) {
      group <- "WT"
    } else if (grepl("KO", folder)) {
      group <- "KO" # 将 Y90C-KO 统一标记为 KO，或者你可以保留 "Y90C_KO"
    } else {
      group <- "Unknown"
    }
    
    # 提取重复编号 (Replicate)
    # 取字符串最后一位作为重复号
    rep_id <- substr(folder, nchar(folder), nchar(folder))
    
    # 写入 Metadata
    sc_obj$Orig_Folder <- folder       # 原始文件夹名
    sc_obj$Group       <- group        # WT vs KO
    sc_obj$Replicate   <- rep_id       # 1 vs 2
    sc_obj$SampleID    <- paste0(group, "_Rep", rep_id) # 比如 WT_Rep1 (用于后续去批次)
    
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-") # 小鼠是 mt-, 人是 MT-
    
    sc_list[[folder]] <- sc_obj
    print(paste("   ✅ 成功入库: Group=", group, "| Rep=", rep_id, "| Cells:", ncol(sc_obj)))
    
  }, error = function(e) {
    message(paste("   ❌ 出错:", folder, e$message))
  })
}

# 合并所有样本
if (length(sc_list) > 0) {
  print("正在合并样本...")
  # add.cell.ids 加上前缀防止细胞条码冲突
  sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = names(sc_list))
} else {
  stop("❌ 未读取到数据")
}
sc_combined <- JoinLayers(sc_combined)
# ------------------------------------------------------------------------------
# 3. 统一质控 (QC)
# ------------------------------------------------------------------------------
print("🚀 步骤2/6: 正在进行质控过滤...")
# 根据你的实际数据情况调整阈值
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# ------------------------------------------------------------------------------
# 4. 去除双细胞 (针对多样本优化版)
# ------------------------------------------------------------------------------
print("🚀 步骤3/6: 识别并去除双细胞 (scDblFinder)...")

# scDblFinder 最佳实践：告诉它哪些细胞属于同一个样本 (samples 参数)
# 这样它会在每个样本内部独立寻找双胞，而不是混在一起找
sce <- as.SingleCellExperiment(sc_combined)
sce <- scDblFinder(sce, samples = "Orig_Folder") # 使用原始文件夹名区分样本

# 将结果导回 Seurat
sc_combined$scDblFinder_class <- sce$scDblFinder.class
n_dbl <- sum(sc_combined$scDblFinder_class == "doublet")
print(paste("   [Result] 总计发现双细胞:", n_dbl, "个"))

# 剔除
sc_combined <- subset(sc_combined, subset = scDblFinder_class == "singlet")

# ------------------------------------------------------------------------------
# 5. 降维与整合 (Integration) - 处理生物学重复的核心
# ------------------------------------------------------------------------------
print("🚀 步骤4/6: 标准化与生物学重复整合 (Harmony)...")

obj <- sc_combined
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)

# --- 关键步骤：使用 Harmony 去除批次效应 (整合 WT_1, WT_2 等) ---
# group.by.vars = "Orig_Folder" 表示我们要消除不同文库(Folder)之间的技术差异
# 这样 WT1 和 WT2 会融合，KO1 和 KO2 会融合
tryCatch({
  obj <- RunHarmony(obj, group.by.vars = "Orig_Folder")
  reduction_to_use <- "harmony"
  print("   ✅ Harmony 整合完成")
}, error = function(e) {
  message("   ⚠️ Harmony 运行失败或未安装，回退到标准 PCA (可能存在批次效应)")
  reduction_to_use <- "pca"
})

# 基于整合后的数据跑 UMAP 和聚类
obj <- RunUMAP(obj, reduction = reduction_to_use, dims = 1:20)
obj <- FindNeighbors(obj, reduction = reduction_to_use, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)

# ------------------------------------------------------------------------------
# 6. ScType 自动注释
# ------------------------------------------------------------------------------
print("🚀 步骤5/6: 运行 ScType 细胞注释...")
# 确保数据库路径正确

db_file_path <- file.path(data_dir, "ScTypeDB_full.xlsx") 

# 如果没有数据库文件，为了防止报错，这里做一个简单的跳过处理
if (file.exists(db_file_path)) {
  gs_list_immune <- gene_sets_prepare(db_file_path, "Brain") # 或者 "All"
  
  # 简化的评分流程
  es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(obj, layer="scale.data")), scaled = TRUE, 
                         gs = gs_list_immune$gs_positive, gs2 = gs_list_immune$gs_negative)
  
  # 映射 Cluster 到 Type
  cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    cells_in_cluster <- rownames(obj@meta.data[obj@meta.data$seurat_clusters == cl, ])
    es.max_subset <- es.max[ , cells_in_cluster, drop = FALSE]
    es.max.cl = sort(rowSums(es.max_subset), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
  }))
  
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # 赋值
  obj@meta.data$cell_type <- ""
  for(j in unique(sctype_scores$cluster)){
    cl_type <- sctype_scores[sctype_scores$cluster==j, "type"]
    obj@meta.data$cell_type[obj@meta.data$seurat_clusters == j] <- as.character(cl_type)
  }
} else {
  message("⚠️ 未找到 ScType 数据库，跳过注释，使用 Cluster ID 绘图")
  obj$cell_type <- obj$seurat_clusters
}
#大脑皮层细胞不应含有Tanycytes，经过鉴定为Mature Neurons
obj[["cell_type"]][obj[["cell_type"]] == "Tanycytes"] <- "Mature neurons"

# ==============================================================================
# 7. 结果可视化与输出 (无清洗版：直接使用原始 ScType 注释 + PNG输出)
# ==============================================================================

print("🚀 步骤6/6: 正在使用原始标签生成 PNG 图片...")

# ------------------------------------------------------------------------------
# A. 设置绘图分组
# ------------------------------------------------------------------------------
# 不进行任何清洗，直接指定使用 obj@meta.data 中的 "cell_type" 列
# 这列包含了 ScType 算出来的原始结果
plot_group <- "cell_type"

print("当前使用的标签列: cell_type (原始未清洗)")

# ------------------------------------------------------------------------------
# B. 准备绘图 (图例移至右侧，防止长标签遮挡)
# ------------------------------------------------------------------------------
# 创建一个新的文件夹存放结果，避免覆盖
plot_dir <- file.path(data_dir, "Results_Plots_Raw") 
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# 通用主题设置
my_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  legend.position = "right",           # 图例放右边
  legend.text = element_text(size = 10), # 图例文字大小
  legend.title = element_blank()       # 去掉图例标题
) 

# 强制图例显示为 1 列
my_guide <- guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

# 1. 总图
p_total <- DimPlot(obj, reduction = "umap", group.by = plot_group, 
                   label = TRUE, label.size = 4, repel = TRUE) + 
  ggtitle(paste0("Total Integrated (Cells: ", ncol(obj), ")")) +
  my_theme + my_guide

# 2. WT 独立图
obj_wt <- subset(obj, subset = Group == "WT")
p_wt <- DimPlot(obj_wt, reduction = "umap", group.by = plot_group, 
                label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("WT Group") +
  my_theme + my_guide

# 3. KO 独立图
obj_ko <- subset(obj, subset = Group == "KO")
p_ko <- DimPlot(obj_ko, reduction = "umap", group.by = plot_group, 
                label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("KO Group") +
  my_theme + my_guide

# 4. 对比图
p_split <- DimPlot(obj, reduction = "umap", group.by = plot_group, split.by = "Group",
                   label = TRUE, label.size = 3, repel = TRUE, ncol = 2) +
  ggtitle("Condition Comparison: WT vs KO") +
  theme(legend.position = "right") + my_guide

# ------------------------------------------------------------------------------
# C. 保存为 PNG 格式
# ------------------------------------------------------------------------------
print(paste("正在保存图片至:", plot_dir))

# 由于原始标签可能很长，这里把宽度(width)设得稍微大一点(14英寸)，防止图例被切掉
ggsave(file.path(plot_dir, "01_UMAP_Total_Raw.png"), plot = p_total, 
       width = 14, height = 9, dpi = 300, bg = "white")

ggsave(file.path(plot_dir, "02_UMAP_WT_Raw.png"), plot = p_wt, 
       width = 14, height = 9, dpi = 300, bg = "white")

ggsave(file.path(plot_dir, "03_UMAP_KO_Raw.png"), plot = p_ko, 
       width = 14, height = 9, dpi = 300, bg = "white")

ggsave(file.path(plot_dir, "04_UMAP_Split_Raw.png"), plot = p_split, 
       width = 20, height = 8, dpi = 300, bg = "white") # 对比图更宽一些

print("✅ 图片生成完毕！请查看 Results_Plots_Raw 文件夹。")
# ==============================================================================
# 8. 细胞比例分析 (带百分比标签)
# ==============================================================================
print("🚀 步骤7/8: 正在进行细胞比例统计分析...")

# 创建统计结果存放文件夹
stats_dir <- file.path(data_dir, "Results_Stats")
if (!dir.exists(stats_dir)) dir.create(stats_dir)

# 8.1 计算比例数据
# ------------------------------------------------------------------------------
# 提取元数据
meta_df <- obj@meta.data

# 使用 dplyr 进行分组统计
prop_data <- meta_df %>%
  group_by(Group, cell_type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(
    Total = sum(Count),
    Percent = Count / Total,
    # 准备标签文本：只有比例大于 3% 才显示数值，防止重叠堆积
    Label = ifelse(Percent > 0.03, paste0(round(Percent * 100, 1), "%"), "") 
  )

# 8.2 绘制堆叠柱状图
# ------------------------------------------------------------------------------
p_barplot <- ggplot(prop_data, aes(x = Group, y = Percent, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  # 添加文字标签 (居中显示)
  geom_text(aes(label = Label), 
            position = position_fill(vjust = 0.5), 
            size = 3.5, color = "black") +
  scale_y_continuous(labels = scales::percent) + # Y轴显示为百分比
  labs(x = "Group", y = "Cell Proportion (%)", title = "Cell Type Proportion: WT vs KO") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_blank()
  )

# 8.3 保存结果
# ------------------------------------------------------------------------------
# 保存图片
ggsave(file.path(stats_dir, "05_Cell_Proportion_Barplot.png"), 
       plot = p_barplot, width = 8, height = 6, dpi = 300)

# 保存统计表格 (CSV)
write.csv(prop_data, file.path(stats_dir, "05_Cell_Proportion_Table.csv"), row.names = FALSE)

print("✅ 细胞比例图已生成！")
# ==============================================================================
# 8.2 [可选] 补充分析：按样本查看重复一致性 (Per Sample Check)
# ==============================================================================
print("🚀 [可选步骤] 正在生成按样本(Sample)拆分的比例图，用于检查重复一致性...")

# 1. 计算按样本的比例数据
# ------------------------------------------------------------------------------
# 注意：group_by 中增加了 SampleID
prop_data_sample <- meta_df %>%
  group_by(SampleID, Group, cell_type) %>%  # 保留 Group 是为了画图时能分栏
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(SampleID) %>%
  mutate(
    Total = sum(Count),
    Percent = Count / Total,
    # 标签逻辑：比例 > 3% 显示数值，否则留空
    Label = ifelse(Percent > 0.03, paste0(round(Percent * 100, 1), "%"), "")
  )

# 2. 绘制分面柱状图 (Facet Plot)
# ------------------------------------------------------------------------------
p_barplot_sample <- ggplot(prop_data_sample, aes(x = SampleID, y = Percent, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_text(aes(label = Label), 
            position = position_fill(vjust = 0.5), 
            size = 3, color = "black") +
  scale_y_continuous(labels = scales::percent) +
  # 关键步骤：使用 facet_grid 按 Group 分栏，让 WT 的两个样本挨在一起，KO 的挨在一起
  facet_grid(~Group, scales = "free_x", space = "free") + 
  labs(x = "Sample ID", y = "Cell Proportion (%)", 
       title = "Cell Type Proportion: Check Replicates Consistency") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # X轴文字斜着放，防止重叠
    legend.position = "right",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#EFEFEF", color = NA), # 分面标题背景色
    strip.text = element_text(face = "bold", size = 12)
  )

# 3. 保存结果
# ------------------------------------------------------------------------------
# 保存图片
ggsave(file.path(stats_dir, "05_Cell_Proportion_Barplot_PerSample.png"), 
       plot = p_barplot_sample, width = 10, height = 6, dpi = 300)

# 保存表格
write.csv(prop_data_sample, file.path(stats_dir, "05_Cell_Proportion_Table_PerSample.csv"), row.names = FALSE)

print("✅ 按样本检查图已生成！请查看: 05_Cell_Proportion_Barplot_PerSample.png")
print("   -> 观察提示：请检查同一组内(如 WT1 和 WT2) 的柱子颜色分布是否高度相似。")

# ==============================================================================
# 9. 差异分析 (KO vs WT, 循环所有细胞类型)
# ==============================================================================
print("🚀 步骤8/8: 正在进行差异基因分析 (KO vs WT)...")

# 9.1 准备工作
# ------------------------------------------------------------------------------
# 切换默认 Idents 为 cell_type，方便后续 subset
Idents(obj) <- "cell_type"

# 获取所有存在的细胞类型
all_cell_types <- unique(obj$cell_type)
deg_list <- list() # 用于存储所有类型的差异基因结果

# 9.2 循环计算
# ------------------------------------------------------------------------------
for (ctype in all_cell_types) {
  print(paste0("  -> 正在分析细胞类型: ", ctype))
  
  # 提取当前细胞类型的子集
  sub_obj <- subset(obj, idents = ctype)
  
  # 切换 Ident 为 Group 以便比较 WT 和 KO
  Idents(sub_obj) <- "Group"
  
  # 检查：必须 WT 和 KO 两个组里都有细胞才能对比
  # table(sub_obj$Group) 会返回各组细胞数
  group_counts <- table(sub_obj$Group)
  
  if (sum(names(group_counts) %in% c("WT", "KO")) == 2 && 
      min(group_counts) >= 3) { # 每个组至少有3个细胞才分析
    
    tryCatch({
      # FindMarkers 计算差异
      # ident.1 = "KO" (实验组), ident.2 = "WT" (基准组)
      # 结果解读：avg_log2FC > 0 代表在 KO 中上调
      markers <- FindMarkers(sub_obj, ident.1 = "KO", ident.2 = "WT", 
                             logfc.threshold = 0.25, # 至少有 0.25 logFC 差异
                             min.pct = 0.1,          # 至少在 10% 细胞中表达
                             only.pos = FALSE)       # 上下调都要
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cell_type <- ctype
        markers$comparison <- "KO_vs_WT"
        
        # 整理列顺序
        markers <- markers %>% select(cell_type, gene, avg_log2FC, p_val_adj, p_val, everything())
        
        # 存入列表
        deg_list[[ctype]] <- markers
      } else {
        message(paste("     [提示] ", ctype, " 没有找到显著差异基因。"))
      }
      
    }, error = function(e) {
      message(paste("     [错误] 分析 ", ctype, " 时出错: ", e$message))
    })
    
  } else {
    message(paste("     [跳过] ", ctype, " 细胞数量不足或缺失某一组 (WT/KO)。"))
  }
}

# 9.3 导出结果
# ------------------------------------------------------------------------------
if (length(deg_list) > 0) {
  # 合并为一个大表格
  all_degs <- do.call(rbind, deg_list)
  
  # 1. 保存所有差异基因到一个 CSV
  write.csv(all_degs, file.path(stats_dir, "06_All_DEGs_KO_vs_WT_Summary.csv"), row.names = FALSE)
  
  # 2. 保存到 Excel (每个细胞类型一个 Sheet，方便查看) - 推荐使用
  wb <- createWorkbook()
  
  # 添加一个总表 Sheet
  addWorksheet(wb, "All_Combined")
  writeData(wb, "All_Combined", all_degs)
  
  # 为每个细胞类型添加单独的 Sheet
  for (ctype_name in names(deg_list)) {
    # Excel Sheet 名字不能太长或包含特殊字符，简单处理一下
    clean_name <- substr(gsub("[^[:alnum:]]", "_", ctype_name), 1, 30)
    addWorksheet(wb, clean_name)
    writeData(wb, clean_name, deg_list[[ctype_name]])
  }
  
  saveWorkbook(wb, file.path(stats_dir, "06_DEGs_KO_vs_WT_Full_Report.xlsx"), overwrite = TRUE)
  
  print(paste("✅ 差异分析完成！结果已保存至:", stats_dir))
  print("   - 06_DEGs_KO_vs_WT_Full_Report.xlsx (推荐查看)")
  
} else {
  print("⚠️ 未找到任何显著差异基因，请检查数据或分组。")
}

print("🎉 所有分析流程结束！")
# # ==============================================================================
# # 7. 结果可视化与输出 (scCustomize 修正版：修复参数报错)
# # ==============================================================================
# 
# # --- 0. 加载必要的包 ---
# if (!require("scCustomize", quietly = TRUE)) {
#   if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#   devtools::install_github("samuel-marsh/scCustomize")
# }
# library(scCustomize)
# library(ggplot2)
# library(scales)
# 
# print("🚀 步骤6/6: 正在使用 scCustomize 生成发表级美图...")
# 
# # ------------------------------------------------------------------------------
# # A. 设置绘图分组与构建稳健的颜色映射
# # ------------------------------------------------------------------------------
# plot_group <- "cell_type" 
# 
# # 1. 获取所有唯一的细胞类型
# unique_types <- sort(unique(as.character(obj@meta.data[[plot_group]])))
# n_types <- length(unique_types)
# 
# print(paste("检测到细胞类型数量:", n_types))
# 
# # 2. 定义柔和的发表级色盘
# my_palette <- c(
#   "#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD", "#802268",
#   "#6BD76B", "#D595A7", "#924822", "#837B8D", "#C75127", "#D58F5C", "#7A65A5", "#E4AF69",
#   "#3B1B53", "#CDDEB7", "#612A79", "#AE1F63", "#E7C453", "#5A655E", "#CC9900", "#99CC00",
#   "#33CC00", "#00CC33", "#00CC99", "#0099CC", "#0033CC", "#3300CC", "#9900CC", "#CC0099",
#   "#CC0033", "#FF3300", "#FF9900", "#FFFF00", "#99FF00", "#33FF00", "#00FF33", "#00FF99",
#   "#0099FF", "#0033FF", "#3300FF", "#9900FF", "#CC00FF", "#FF00CC", "#FF0033", "#FF3333"
# )
# 
# # 3. 截取并绑定名字
# if(n_types > length(my_palette)){
#   final_colors <- scales::hue_pal()(n_types)
# } else {
#   final_colors <- my_palette[1:n_types]
# }
# names(final_colors) <- unique_types 
# 
# # ------------------------------------------------------------------------------
# # B. 定义增强版箭头主题 (图例位置在这里控制)
# # ------------------------------------------------------------------------------
# arrow_theme <- theme(
#   axis.line = element_line(arrow = arrow(length = unit(0.25, "cm"), type = "closed"), size = 1), 
#   axis.title = element_text(size = 14, face = "bold", hjust = 0.05), 
#   plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), 
#   legend.text = element_text(size = 12),
#   legend.position = "right" # 【关键】图例位置必须写在 theme 里
# )
# 
# # ------------------------------------------------------------------------------
# # C. 绘图与保存
# # ------------------------------------------------------------------------------
# plot_dir <- file.path(data_dir, "Results_Plots_scCustomize") 
# if (!dir.exists(plot_dir)) dir.create(plot_dir)
# 
# print(paste("正在保存图片至:", plot_dir))
# 
# # --- 1. Total 图 ---
# print("正在绘制: Total Integrated...")
# p_total <- DimPlot_scCustom(
#   seurat_object = obj, 
#   group.by = plot_group, 
#   colors_use = final_colors,  
#   figure_plot = TRUE,         
#   label = FALSE,              
#   pt.size = 0.8               
#   # 【修复】这里删除了 legend.position 参数
# ) + arrow_theme + ggtitle(paste0("Total (Cells: ", ncol(obj), ")"))
# 
# ggsave(file.path(plot_dir, "01_UMAP_Total_scCustom.png"), p_total, width = 14, height = 12, dpi = 300)
# 
# # --- 2. WT 独立图 ---
# print("正在绘制: WT Group...")
# obj_wt <- subset(obj, subset = Group == "WT")
# p_wt <- DimPlot_scCustom(
#   seurat_object = obj_wt, 
#   group.by = plot_group, 
#   colors_use = final_colors,  
#   figure_plot = TRUE,
#   label = FALSE,
#   pt.size = 0.8
# ) + arrow_theme + ggtitle("WT Group")
# 
# ggsave(file.path(plot_dir, "02_UMAP_WT_scCustom.png"), p_wt, width = 14, height = 12, dpi = 300)
# 
# # --- 3. KO 独立图 ---
# print("正在绘制: KO Group...")
# obj_ko <- subset(obj, subset = Group == "KO")
# p_ko <- DimPlot_scCustom(
#   seurat_object = obj_ko, 
#   group.by = plot_group, 
#   colors_use = final_colors,
#   figure_plot = TRUE,
#   label = FALSE,
#   pt.size = 0.8
# ) + arrow_theme + ggtitle("KO Group")
# 
# ggsave(file.path(plot_dir, "03_UMAP_KO_scCustom.png"), p_ko, width = 14, height = 12, dpi = 300)
# 
# # --- 4. 对比图 (Split View) ---
# print("正在绘制: Split Comparison...")
# p_split <- DimPlot_scCustom(
#   seurat_object = obj, 
#   group.by = plot_group, 
#   split.by = "Group",         
#   colors_use = final_colors,
#   figure_plot = TRUE,
#   label = FALSE,
#   pt.size = 0.8,
#   num_columns = 2            
# ) + arrow_theme + ggtitle("Condition Comparison: WT vs KO")
# 
# ggsave(file.path(plot_dir, "04_UMAP_Split_scCustom.png"), p_split, width = 16, height = 8, dpi = 300)
# 
# print("✅ 修复完成！图片已成功生成。")
