# ==============================================================================
# Aorta 单组织分析流程 (Seurat v5) - 包含：加载/QC/去双胞/注释/绘图/差异分析
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(openxlsx) 
  library(Matrix)   
  library(HGNChelper)
  library(scales)
  library(scDblFinder)
  library(SingleCellExperiment)
})

# 加载 ScType 核心函数
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# ------------------------------------------------------------------------------
# 1. 设置路径与样本映射
# ------------------------------------------------------------------------------
root_dir <- "/mnt/disk1/qiuzerui/coldmouse"
out_dir  <- file.path(root_dir, "results")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 本地数据库路径
db_file_path <- file.path(root_dir, "ScTypeDB_full.xlsx") 

# 定义样本与分组
rt_samples <- c("SRR35688257", "SRR35688258", "SRR35688259")
sample_ids <- c(rt_samples, "SRR35688260", "SRR35688261", "SRR35688262")

# ------------------------------------------------------------------------------
# 2. 读取数据并添加 Metadata
# ------------------------------------------------------------------------------
print("🚀 步骤1: 加载 6 个样本数据...")
sc_list <- list()

for (sample in sample_ids) {
  matrix_path <- file.path(root_dir, paste0("Output_", sample), "outs/filtered_feature_bc_matrix")
  if(!dir.exists(matrix_path)) next
  
  counts <- Read10X(data.dir = matrix_path)
  sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
  
  # 注入分组信息
  sc_obj$Group <- ifelse(sample %in% rt_samples, "RT", "LT")
  sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
  sc_list[[sample]] <- sc_obj
}

# 合并对象
sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = sample_ids)

# ------------------------------------------------------------------------------
# 3. 质控与图层合并 (Seurat v5 关键步)
# ------------------------------------------------------------------------------
print("🚀 步骤2: 质控与图层合并...")

# 过滤线粒体 > 15% 的细胞
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# 【V5 核心】合并图层以支持后续全局分析
sc_combined <- JoinLayers(sc_combined) 

# 移除核糖体基因
non_ribo_genes <- setdiff(rownames(sc_combined), grep("^Rp[sl]", rownames(sc_combined), value = T, ignore.case = T))
sc_combined <- subset(sc_combined, features = non_ribo_genes)

# ------------------------------------------------------------------------------
# 4. 去双胞与标准流水线
# ------------------------------------------------------------------------------
print("🚀 步骤3: 运行去双胞与 UMAP 降维...")

# 基于 Join 后的矩阵运行 scDblFinder
sce <- as.SingleCellExperiment(sc_combined)
sce <- scDblFinder(sce, samples = "orig.ident") 
sc_combined$scDblFinder_class <- sce$scDblFinder.class
sc_combined <- subset(sc_combined, subset = scDblFinder_class == "singlet")

# 标准流程
sc_combined <- NormalizeData(sc_combined) %>% 
               FindVariableFeatures(nfeatures = 2000) %>% 
               ScaleData() %>% 
               RunPCA(verbose = FALSE) %>% 
               RunUMAP(dims = 1:20) %>% 
               FindNeighbors(dims = 1:20) %>% 
               FindClusters(resolution = 0.5)

# ------------------------------------------------------------------------------
# 5. ScType 自动注释 (修正版：适配小鼠 & 主动脉复杂组织)
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在使用本地数据库进行细胞注释...")

tryCatch({
  # 【关键修改 1】扩充组织类型
  # 主动脉不仅有免疫细胞，还有平滑肌(Muscle)、内皮(Endothelial)等
  # 注意：ScTypeDB 中通常用 "Immune system", "Muscle", "Heart" 等大类
  target_tissues <- c("Immune system") 
  
  gs_list <- gene_sets_prepare(db_file_path, target_tissues)
  
  # 提取矩阵
  sc_data <- as.matrix(LayerData(sc_combined, layer = "scale.data")) 
  
  # 【关键修改 2】强制转换为大写以适配小鼠数据
  # 小鼠基因是 Title Case (Actb)，数据库是 Upper Case (ACTB)
  rownames(sc_data) <- toupper(rownames(sc_data))
  print("   ℹ️ 已将基因名转换为大写以适配 ScType 数据库")

  # 计算打分
  es.max <- sctype_score(scRNAseqData = sc_data, scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # 按 Cluster 汇总结果
  cL_resutls <- do.call("rbind", lapply(unique(sc_combined$seurat_clusters), function(cl){
    cells_in_cluster <- WhichCells(sc_combined, idents = cl)
    valid_cells <- intersect(colnames(es.max), cells_in_cluster)
    es.max_subset <- es.max[ , valid_cells, drop = FALSE]
    
    if(ncol(es.max_subset) > 0) {
      best_type <- names(sort(rowSums(es.max_subset), decreasing = TRUE))[1]
    } else {
      best_type <- "Unknown"
    }
    data.frame(cluster = cl, type = best_type)
  }))
  
  # 建立 Cluster -> Type 的映射
  cluster_to_type <- setNames(cL_resutls$type, cL_resutls$cluster)
  
  # 赋值给 metadata (使用 unname 强制按顺序赋值)
  current_clusters <- as.character(sc_combined$seurat_clusters)
  new_types <- cluster_to_type[current_clusters]
  sc_combined$cell_type <- unname(new_types)
  
  print("✅ 注释成功！")
  print(table(sc_combined$cell_type)) # 打印一下看看有没有平滑肌细胞
  
}, error = function(e) {
  message("❌ 注释失败，回退到数字编号。详细错误: ", e$message)
  sc_combined$cell_type <<- unname(as.character(sc_combined$seurat_clusters))
})

# 再次检查
if(!"cell_type" %in% colnames(sc_combined@meta.data)){
  sc_combined$cell_type <- unname(as.character(sc_combined$seurat_clusters))
}

# ------------------------------------------------------------------------------
# 6. 绘图 (重构版：独立输出三张图)
# ------------------------------------------------------------------------------
print("🚀 步骤5: 开始绘图 (独立文件输出)...")

# 定义当前组织名称
tissue_name <- "Aorta"

# 安全检查
if(!"cell_type" %in% colnames(sc_combined@meta.data)) {
  sc_combined$cell_type <- as.character(sc_combined$seurat_clusters)
}

tryCatch({
  # --- 1. 锁定因子水平 (核心：保证三张图颜色一致) ---
  all_cell_types <- sort(unique(sc_combined$cell_type))
  sc_combined$cell_type <- factor(sc_combined$cell_type, levels = all_cell_types)
  
  # --- 2. 构建颜色字典 (Named Vector) ---
  my_colors <- hue_pal()(length(all_cell_types))
  names(my_colors) <- all_cell_types
  
  # --- 3. 绘制并保存：Total 图 ---
  print("   📸 正在绘制并保存 Total 层...")
  p_total <- DimPlot(sc_combined, reduction = "umap", group.by = "cell_type", 
                     cols = my_colors, 
                     label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "- Total")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 保存 Total
  ggsave(filename = paste0(tissue_name, "_Total.png"), 
         plot = p_total, width = 10, height = 8, path = out_dir)
  
  # --- 4. 绘制并保存：LT 分组 ---
  print("   📸 正在绘制并保存 LT 分组...")
  # 注意：这里去掉了 NoLegend()，因为独立图片必须要有图例
  p_lt <- DimPlot(subset(sc_combined, subset = Group == "LT"), 
                  reduction = "umap", group.by = "cell_type", 
                  cols = my_colors, 
                  label = FALSE) +  # 分组图通常不加标签，防遮挡，通过图例看
    ggtitle("LT (Cold)") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 保存 LT
  ggsave(filename = paste0(tissue_name, "_LT.png"), 
         plot = p_lt, width = 10, height = 8, path = out_dir)
  
  # --- 5. 绘制并保存：RT 分组 ---
  print("   📸 正在绘制并保存 RT 分组...")
  p_rt <- DimPlot(subset(sc_combined, subset = Group == "RT"), 
                  reduction = "umap", group.by = "cell_type", 
                  cols = my_colors, 
                  label = FALSE) + 
    ggtitle("RT (Room Temp)") + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 保存 RT
  ggsave(filename = paste0(tissue_name, "_RT.png"), 
         plot = p_rt, width = 10, height = 8, path = out_dir)

  print("✅ 所有图片已独立保存成功！")

}, error = function(e) {
  print(paste("❌ 绘图阶段发生错误:", e$message))
})
# 7. 差异表达分析 (LT vs RT)
# ------------------------------------------------------------------------------
print("🚀 步骤6: 正在执行 LT vs RT 差异分析...")

de_dir <- file.path(out_dir, "DE_Results")
if(!dir.exists(de_dir)) dir.create(de_dir)

for (ctype in unique(sc_combined$cell_type)) {
  tryCatch({
    sub_obj <- subset(sc_combined, subset = cell_type == ctype)
    if(sum(sub_obj$Group == "LT") >= 3 && sum(sub_obj$Group == "RT") >= 3) {
      Idents(sub_obj) <- "Group"
      markers <- FindMarkers(sub_obj, ident.1 = "LT", ident.2 = "RT", logfc.threshold = 0.25)
      markers$gene <- rownames(markers)
      write.xlsx(markers, file = file.path(de_dir, paste0("DE_", gsub("/", "_", ctype), ".xlsx")))
    }
  }, error = function(e) next)
}

saveRDS(sc_combined, file = file.path(out_dir, "Aorta_Final_Object.rds"))
print("🎉 全部任务顺利完成！")
# ==============================================================================
# 8. 细胞类型比例分析 (堆积柱状图)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(scales) # 用于百分比格式化

print("🚀 步骤8: 正在绘制细胞比例堆积柱状图 (带标签)...")

# 8.0 准备工作：设置分组因子顺序
# 确保 RT 在前，LT 在后 (或者根据你的需求调整顺序)
if("Group" %in% colnames(sc_combined@meta.data)){
  sc_combined$Group <- factor(sc_combined$Group, levels = c("RT", "LT"))
}

# 8.1 构建绘图数据
# 统计每个 Group 中每种 cell_type 的细胞数量
cell_stats <- sc_combined@meta.data %>%
  group_by(Group, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(percentage = count / sum(count)) # 计算小数比例

# --- 准备标签文本 ---
# 使用 scales::percent 将小数转换为百分比字符串，保留1位小数
cell_stats$label_text <- scales::percent(cell_stats$percentage, accuracy = 0.1)

# 确保 cell_type 因子顺序与之前一致（为了颜色对应）
# 这里复用 Step 6 生成的 my_colors 的名字作为 level 顺序
if(exists("my_colors")){
  cell_stats$cell_type <- factor(cell_stats$cell_type, levels = names(my_colors))
}

# 8.2 绘图
p_bar <- ggplot(cell_stats, aes(x = Group, y = percentage, fill = cell_type)) +
  # 绘制柱状图层
  geom_bar(stat = "identity", position = "fill", width = 0.7) + 
  
  # --- 添加文字标签图层 ---
  # 逻辑：只显示比例 > 1% (0.01) 的标签，防止字重叠
  geom_text(data = subset(cell_stats, percentage > 0.01), 
            aes(label = label_text),
            position = position_fill(vjust = 0.5),
            size = 3,          # 字体大小
            color = "black") + # 字体颜色
  
  # 设置Y轴格式
  scale_y_continuous(labels = scales::percent) + 
  
  # --- 关键：使用与 UMAP 一致的颜色 ---
  # 如果 Step 6 定义了 my_colors，则使用它；否则使用默认色
  (if(exists("my_colors")) scale_fill_manual(values = my_colors) else scale_fill_hue()) +
  
  labs(title = "Aorta Cell Type Proportion by Group",
       x = "Condition", 
       y = "Percentage of Cells",
       fill = "Cell Type") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

# 8.3 保存
filename_bar <- "Aorta_CellType_Proportion_Labeled.png"
ggsave(filename = filename_bar, plot = p_bar, width = 8, height = 6, path = out_dir)

print(paste("   ✅ 带标签的比例图已保存:", filename_bar))
