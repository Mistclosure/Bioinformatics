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
# 5. ScType 自动注释 (修正版)
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在使用本地数据库进行细胞注释...")

tryCatch({
  # 准备基因集
  gs_list <- gene_sets_prepare(db_file_path, "Immune system")
  
  # 【修复1】确保输入是普通矩阵，Seurat V5 的 LayerData 有时是稀疏矩阵
  # 注意：如果内存不足，可以分批处理，但 ScType 需要密集矩阵计算
  sc_data <- as.matrix(LayerData(sc_combined, layer = "scale.data")) 
  
  # 计算打分
  es.max <- sctype_score(scRNAseqData = sc_data, scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # 按 Cluster 汇总结果
  cL_resutls <- do.call("rbind", lapply(unique(sc_combined$seurat_clusters), function(cl){
    cells_in_cluster <- WhichCells(sc_combined, idents = cl)
    # 确保只取存在的细胞
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
  
  # 【关键修复2】使用 unname() 去除名称，防止 "No cell overlap" 错误
  # 将 Cluster ID 映射为 Cell Type
  current_clusters <- as.character(sc_combined$seurat_clusters)
  new_types <- cluster_to_type[current_clusters]
  
  # 赋值给 metadata (使用 unname 强制按顺序赋值)
  sc_combined$cell_type <- unname(new_types)
  
  print("✅ 注释成功！")
  
}, error = function(e) {
  message("❌ 注释失败，回退到数字编号。详细错误: ", e$message)
  # 【关键修复3】错误处理中也要加 unname()
  sc_combined$cell_type <<- unname(as.character(sc_combined$seurat_clusters))
})

# 检查 cell_type 是否成功创建，防止绘图报错
if(!"cell_type" %in% colnames(sc_combined@meta.data)){
  message("⚠️ 警告: cell_type 未创建，强制使用聚类编号")
  sc_combined$cell_type <- unname(as.character(sc_combined$seurat_clusters))
}

# ------------------------------------------------------------------------------
# 6. 绘图结果导出 (修正版)
# ------------------------------------------------------------------------------
print("🚀 步骤5: 生成结果图表...")

# 确保 cell_type 是因子或字符
sc_combined$cell_type <- as.factor(sc_combined$cell_type)

# 锁定颜色映射
all_cell_types <- sort(unique(sc_combined$cell_type))
my_colors <- hue_pal()(length(all_cell_types))
names(my_colors) <- all_cell_types

p1 <- DimPlot(sc_combined, reduction = "umap", group.by = "cell_type", cols = my_colors, label = TRUE) + 
      ggtitle("Aorta Annotation") + NoLegend()

p2 <- DimPlot(sc_combined, reduction = "umap", group.by = "Group", cols = c("RT" = "#A6CEE3", "LT" = "#1F78B4")) + 
      ggtitle("RT vs LT Distribution")

ggsave("Aorta_Annotation_Final.png", plot = p1 + p2, path = out_dir, width = 14, height = 6)

# ------------------------------------------------------------------------------
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
