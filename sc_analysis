# ==============================================================================
# 完整流程代码：加载 -> 质控 -> 去双胞 -> 注释 -> 绘图(修正版)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx) 
library(Matrix)   
library(HGNChelper)
library(scales) # 用于生成默认色盘 (绘图必需)

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
# 1. 设置路径与元数据映射
# ------------------------------------------------------------------------------
# !!! 请务必修改为您电脑上的实际路径 !!!
data_dir <- "D:/qiuzerui/单细胞寒冷处理小鼠" 

sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")
tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")
group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")

# ------------------------------------------------------------------------------
# 2. 读取数据、修复格式并添加 Metadata
# ------------------------------------------------------------------------------
sc_list <- list()
full_folder_names <- list.files(data_dir, pattern = "_matrix")

print("🚀 步骤1/6: 开始读取并修复 Singleron 数据...")

for (sample in sample_ids) {
  # 模糊匹配文件夹名
  folder <- full_folder_names[grep(paste0("^", sample, "_"), full_folder_names)]
  
  if(length(folder) == 0) {
    message(paste("⚠️ 跳过：未找到样本", sample)); next
  }
  
  print(paste("   正在处理:", sample))
  
  tryCatch({
    counts <- Read10X(data.dir = file.path(data_dir, folder[1]))
    
    if (is.list(counts) && !is(counts, "dgCMatrix")) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts$`Gene Expression` 
      } else {
        counts <- counts[[1]]
      }
    }
    
    colnames(counts) <- gsub("_", "-", colnames(counts))
    
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    if (ncol(sc_obj) == 0) { message(paste("   ⚠️ 警告：样本无细胞，跳过。")); next }
    
    # 解析元数据
    prefix <- substr(sample, 1, 1)
    suffix <- substr(sample, 2, 2)
    
    sc_obj$Tissue <- as.character(tissue_map[prefix])
    sc_obj$Group  <- as.character(group_map[suffix])
    
    if(suffix == "2") {
      sc_obj$Condition <- "Cold"
    } else {
      sc_obj$Condition <- "NonCold"
    }
    
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
    
    sc_list[[sample]] <- sc_obj
    print(paste("   ✅ 成功入库: ", sample, "| 细胞数:", ncol(sc_obj)))
    
  }, error = function(e) {
    message(paste("   ❌ 处理样本", sample, "时出错:", e$message))
  })
}

if (length(sc_list) > 0) {
  print("正在合并所有样本...")
  sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = sample_ids)
} else {
  stop("❌ 未读取到任何数据！请检查路径。")
}

# ------------------------------------------------------------------------------
# 3. 统一质控 (QC)
# ------------------------------------------------------------------------------
print("🚀 步骤2/6: 正在进行质控过滤...")
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# ------------------------------------------------------------------------------
# 4. 按组织拆分并独立分析 (含 scDblFinder 去双胞)
# ------------------------------------------------------------------------------
print("🚀 步骤3/6: 按组织拆分并独立聚类 (含 scDblFinder 去双胞)...")

sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

# 定义标准处理函数
run_standard_pipeline <- function(obj, tissue_name) {
  print(paste(">>> 正在处理组织:", tissue_name, "| 初始细胞数:", ncol(obj)))
  
  if(ncol(obj) < 50) { return(NULL) }
  
  # A. 基础预处理 (为了给 scDblFinder 提供聚类信息)
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5) 
  
  # B. 运行 scDblFinder 去除双细胞
  print(paste("   [scDblFinder] 正在识别双细胞..."))
  
  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class
    
    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    print(paste("   [Result] 发现双细胞:", n_dbl, "个 (占比", round(n_dbl/ncol(obj)*100, 2), "%)"))
    
    # --- 剔除双细胞 ---
    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    print(paste("   [Filter] 剔除后剩余细胞:", ncol(obj)))
    
    # C. 剔除后重新进行特征选择和标准化
    print("   [Re-Process] 正在基于纯净细胞重新寻找高变基因...")
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:20)
    
  }, error = function(e) {
    message(paste("   ⚠️ scDblFinder 运行警告:", e$message, "- 将保留所有细胞继续"))
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:20)
  })
  
  return(obj)
}

# 批量运行 pipeline
sc_by_tissue <- lapply(names(sc_by_tissue), function(x) {
  run_standard_pipeline(sc_by_tissue[[x]], x)
})

names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow") 
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]

# ------------------------------------------------------------------------------
# 5. ScType 自动细胞注释
# ------------------------------------------------------------------------------
print("🚀 步骤4/6: 运行 ScType 细胞注释...")

db_file_path <- file.path(data_dir, "ScTypeDB_full.xlsx") 

if (!file.exists(db_file_path)) {
  stop(paste("❌ 未找到数据库文件！请确保 ScTypeDB_full.xlsx 在路径:", data_dir))
}

gs_list_immune <- gene_sets_prepare(db_file_path, "Immune system") 

# 注释函数
run_annotation <- function(obj, gs_list, custom_name = NULL) {
  
  if (!is.null(custom_name)) { obj_name_str <- custom_name } else { obj_name_str <- deparse(substitute(obj)) }
  
  # ScType 评分
  es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(obj, layer="scale.data")), scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # 映射到 Cluster
  cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    cells_in_cluster <- rownames(obj@meta.data[obj@meta.data$seurat_clusters == cl, ])
    es.max_subset <- es.max[ , cells_in_cluster, drop = FALSE] # 关键修复
    es.max.cl = sort(rowSums(es.max_subset), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
  }))
  
  # 写入 Metadata
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  obj@meta.data$cell_type <- ""
  for(j in unique(sctype_scores$cluster)){
    cl_type <- sctype_scores[sctype_scores$cluster==j, "type"]
    obj@meta.data$cell_type[obj@meta.data$seurat_clusters == j] <- as.character(cl_type)
  }
  return(obj)
}

# 批量运行注释
sc_by_tissue_annotated <- lapply(names(sc_by_tissue), function(nm) {
  current_obj <- sc_by_tissue[[nm]]
  new_obj <- run_annotation(current_obj, gs_list_immune, custom_name = nm)
  return(new_obj)
})

names(sc_by_tissue_annotated) <- names(sc_by_tissue)
sc_by_tissue <- sc_by_tissue_annotated

# ------------------------------------------------------------------------------
# 6. 绘图 (最终修正版：单图例 + 颜色锁定)
# ------------------------------------------------------------------------------
print("🚀 步骤5/6: 开始绘图 (已修复重复图例问题)...")

for (tissue_name in names(sc_by_tissue)) {
  
  print(paste("   正在绘制:", tissue_name, "..."))
  obj <- sc_by_tissue[[tissue_name]]
  
  # 1. 锁定因子水平
  all_cell_types <- sort(unique(obj$cell_type))
  obj$cell_type <- factor(obj$cell_type, levels = all_cell_types)
  
  # 2. 构建颜色字典 (Named Vector)
  # 确保子图即使缺失某些细胞类型，颜色也保持一致
  my_colors <- hue_pal()(length(all_cell_types))
  names(my_colors) <- all_cell_types
  
  # 3. 绘制 Total 图 (保留图例)
  p_total <- DimPlot(obj, reduction = "umap", group.by = "cell_type", 
                     cols = my_colors, 
                     label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "- Total")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 4. 绘制各分组子图 (强制隐藏图例 NoLegend)
  p_cold <- DimPlot(subset(obj, subset = Group == "Cold_4C"), 
                    reduction = "umap", group.by = "cell_type", 
                    cols = my_colors, 
                    label = FALSE) + 
    ggtitle("Cold_4C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    NoLegend() # 关键
  
  p_rt <- DimPlot(subset(obj, subset = Group == "RT_25C"), 
                  reduction = "umap", group.by = "cell_type", 
                  cols = my_colors, 
                  label = FALSE) + 
    ggtitle("RT_25C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    NoLegend() # 关键
  
  p_tn <- DimPlot(subset(obj, subset = Group == "TN_30C"), 
                  reduction = "umap", group.by = "cell_type", 
                  cols = my_colors, 
                  label = FALSE) + 
    ggtitle("TN_30C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    NoLegend() # 关键
  
  # 5. 拼图
  # 注意：去掉了 legend.position = "right"，让 patchwork 自动使用 p_total 的图例
  p_final <- (p_total | p_cold) / (p_rt | p_tn) + 
    plot_layout(guides = "collect", axes = "collect") &
    theme(
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "cm")
    )
  
  # 6. 保存
  filename_png <- paste0(tissue_name, "_Annotation_Fixed_Final.png")
  
  tryCatch({
    ggsave(filename = filename_png, plot = p_final, width = 16, height = 12, path = data_dir)
    print(paste("   ✅ [", tissue_name, "] 图片已保存:", filename_png))
  }, error = function(e) {
    print(paste("   ❌ [", tissue_name, "] 保存失败:", e$message))
  })
}

print("🎉 全部任务完成！")
# ==============================================================================
# 7. 差异表达分析 (DE Analysis) - 文件夹分层版
# ==============================================================================
print("🚀 步骤6/6: 开始差异分析 (按文件夹分层存储)...")

library(openxlsx)
library(dplyr)
library(Seurat)

# 1. 定义总输出目录
de_root_dir <- file.path(data_dir, "DE_Results")
if(!dir.exists(de_root_dir)) dir.create(de_root_dir)

# 2. 循环处理每个组织
for (tissue_name in names(sc_by_tissue)) {
  
  print(paste(">>> 正在分析组织:", tissue_name))
  obj <- sc_by_tissue[[tissue_name]]
  DefaultAssay(obj) <- "RNA" # 确保使用 RNA assay
  
  # --- 创建当前组织的两个对比文件夹 ---
  # 路径: DE_Results/Aorta/Cold_vs_RT_25C
  dir_comp1 <- file.path(de_root_dir, tissue_name, "Cold_vs_RT_25C")
  # 路径: DE_Results/Aorta/Cold_vs_TN_30C
  dir_comp2 <- file.path(de_root_dir, tissue_name, "Cold_vs_TN_30C")
  
  if(!dir.exists(dir_comp1)) dir.create(dir_comp1, recursive = TRUE)
  if(!dir.exists(dir_comp2)) dir.create(dir_comp2, recursive = TRUE)
  
  cell_types <- unique(obj$cell_type)
  
  # 3. 循环处理每个细胞类型
  for (ctype in cell_types) {
    # 格式化文件名安全字符
    safe_ctype_name <- gsub("[[:punct:]]", "_", ctype)
    # 文件名前缀: 组织_细胞类型
    file_prefix <- paste0(tissue_name, "_", safe_ctype_name)
    
    print(paste("   --- 正在计算:", ctype, "---"))
    
    tryCatch({
      sub_obj <- subset(obj, subset = cell_type == ctype)
      Idents(sub_obj) <- "Group"
      
      # 统计细胞数
      n_cold <- sum(sub_obj$Group == "Cold_4C")
      n_rt   <- sum(sub_obj$Group == "RT_25C")
      n_tn   <- sum(sub_obj$Group == "TN_30C")
      
      # ==========================================================
      # 对比组 1: Cold_4C vs RT_25C (保存到 dir_comp1)
      # ==========================================================
      if (n_cold >= 3 && n_rt >= 3) {
        markers_rt <- FindMarkers(sub_obj, ident.1 = "Cold_4C", ident.2 = "RT_25C", 
                                  logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
        
        if (nrow(markers_rt) > 0) {
          # 整理表格
          markers_rt$gene <- rownames(markers_rt)
          markers_rt$comparison <- "Cold_vs_RT"
          markers_rt$cell_type <- ctype
          markers_rt <- markers_rt %>% select(gene, everything())
          
          # 保存文件
          fname <- paste0(file_prefix, ".xlsx")
          write.xlsx(markers_rt, file = file.path(dir_comp1, fname))
          print(paste("      ✅ [VS RT] 保存:", fname))
        } else {
          print("      ⚠️ [VS RT] 无差异基因")
        }
      } else {
        print(paste0("      ⏭️ [VS RT] 跳过 (细胞不足: ", n_cold, " vs ", n_rt, ")"))
      }
      
      # ==========================================================
      # 对比组 2: Cold_4C vs TN_30C (保存到 dir_comp2)
      # ==========================================================
      if (n_cold >= 3 && n_tn >= 3) {
        markers_tn <- FindMarkers(sub_obj, ident.1 = "Cold_4C", ident.2 = "TN_30C", 
                                  logfc.threshold = 0.25, min.pct = 0.1, only.pos = FALSE)
        
        if (nrow(markers_tn) > 0) {
          # 整理表格
          markers_tn$gene <- rownames(markers_tn)
          markers_tn$comparison <- "Cold_vs_TN"
          markers_tn$cell_type <- ctype
          markers_tn <- markers_tn %>% select(gene, everything())
          
          # 保存文件
          fname <- paste0(file_prefix, ".xlsx")
          write.xlsx(markers_tn, file = file.path(dir_comp2, fname))
          print(paste("      ✅ [VS TN] 保存:", fname))
        } else {
          print("      ⚠️ [VS TN] 无差异基因")
        }
      } else {
        print(paste0("      ⏭️ [VS TN] 跳过 (细胞不足: ", n_cold, " vs ", n_tn, ")"))
      }
      
    }, error = function(e) {
      print(paste("   ❌ 出错:", ctype, e$message))
    })
  }
}

print("🎉 分析完成！请查看 DE_Results 文件夹内的分类结果。")
