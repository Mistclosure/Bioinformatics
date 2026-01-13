# ==============================================================================
# 完整流程代码：加载 -> 质控 -> 去双胞 -> 注释 -> 绘图(修正版)->差异分析
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
# 3. 统一质控 (QC) 与 移除核糖体基因 (新增)
# ------------------------------------------------------------------------------
print("🚀 步骤2/6: 正在进行质控过滤并移除核糖体基因...")

# A. 细胞水平过滤 (保持原样)
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# B. 计算核糖体基因占比 (可选，仅为了查看过滤前的情况)
# 小鼠核糖体基因通常以 Rp[sl] 开头
sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = "^Rp[sl]")

# C. 基因水平过滤：彻底从矩阵中移除核糖体基因
# 匹配所有以 Rps 或 Rpl 开头的基因 (不区分大小写以防万一)
ribo_genes <- grep("^Rp[sl]", rownames(sc_combined), value = TRUE, ignore.case = TRUE)

# 打印移除信息，方便检查
print(paste("📊 识别到核糖体基因数量:", length(ribo_genes)))
print(paste("🧬 过滤前总基因数:", nrow(sc_combined)))

# 提取非核糖体基因的名称
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]

# 执行过滤：仅保留非核糖体基因
sc_combined <- subset(sc_combined, features = non_ribo_genes)

print(paste("✅ 过滤后剩余基因数:", nrow(sc_combined)))
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
# ==============================================================================
# 5. ScType 自动细胞注释 (修正版：组织特异性 + PBMC 强制移除 Macrophage)
# ==============================================================================
print("🚀 步骤4/6: 运行 ScType 细胞注释 (已应用组织特异性策略)...")

db_file_path <- file.path(data_dir, "ScTypeDB_full.xlsx") 

if (!file.exists(db_file_path)) {
  stop(paste("❌ 未找到数据库文件！请确保 ScTypeDB_full.xlsx 在路径:", data_dir))
}

# 定义注释执行函数 (保持不变)
run_annotation <- function(obj, gs_list, custom_name = NULL) {
  
  if (!is.null(custom_name)) { obj_name_str <- custom_name } else { obj_name_str <- deparse(substitute(obj)) }
  
  # ScType 评分
  es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(obj, layer="scale.data")), scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # 映射到 Cluster
  cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    cells_in_cluster <- rownames(obj@meta.data[obj@meta.data$seurat_clusters == cl, ])
    es.max_subset <- es.max[ , cells_in_cluster, drop = FALSE] 
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

# --- 关键修改：针对不同组织准备不同的 Gene Sets ---
sc_by_tissue_annotated <- lapply(names(sc_by_tissue), function(nm) {
  
  current_obj <- sc_by_tissue[[nm]]
  print(paste(">>> 正在为", nm, "准备 ScType 注释..."))
  
  # 1. 确定组织数据库策略##########针对不同组织修改！！！！！！！！！######
  if (nm == "Aorta") {
    target_tissues <- c("Immune system")
    print("   -> 策略: Aorta (Immune)")
    
  } else if (nm == "BoneMarrow") {
    target_tissues <- c("Immune system")
    print("   -> 策略: BoneMarrow (Immune )")
    
  } else {
    # PBMC 处理
    target_tissues <- c("Immune system")
    print("   -> 策略: PBMC (仅 Immune system)")
  }
  
  # 2. 动态生成 gs_list 并进行特殊过滤
  tryCatch({
    # 加载原始基因集
    gs_list_dynamic <- gene_sets_prepare(db_file_path, target_tissues)
    
    # --- 【新增】针对 PBMC 的特殊过滤 ---
    if (nm == "PBMC") {
      # 查找所有名字里包含 "Macrophage" 的细胞类型 (不区分大小写)
      types_to_remove <- grep("Macrophage", names(gs_list_dynamic$gs_positive), ignore.case = TRUE, value = TRUE)
      
      if (length(types_to_remove) > 0) {
        print(paste("   -> 🛑 [PBMC特异性修正] 正在移除巨噬细胞选项:", paste(types_to_remove, collapse = ", ")))
        
        # 从正向和负向列表中移除这些细胞
        gs_list_dynamic$gs_positive <- gs_list_dynamic$gs_positive[ !names(gs_list_dynamic$gs_positive) %in% types_to_remove ]
        gs_list_dynamic$gs_negative <- gs_list_dynamic$gs_negative[ !names(gs_list_dynamic$gs_negative) %in% types_to_remove ]
      }
    }
    # ------------------------------------
    
    # 执行注释
    new_obj <- run_annotation(current_obj, gs_list_dynamic, custom_name = nm)
    return(new_obj)
    
  }, error = function(e) {
    message(paste("   ❌ 准备数据库或注释时出错:", e$message))
    return(current_obj) 
  })
})

names(sc_by_tissue_annotated) <- names(sc_by_tissue)
sc_by_tissue <- sc_by_tissue_annotated
# ==============================================================================
# 5.5 特殊处理：合并 PBMC 中的单核细胞亚群 (新增)
# ==============================================================================
if ("PBMC" %in% names(sc_by_tissue)) {
  print("🚀 正在执行 PBMC 特殊处理：合并 Monocytes 亚群...")
  
  # 提取 PBMC 对象
  pbmc_obj <- sc_by_tissue[["PBMC"]]
  
  # 记录合并前的类型
  old_types <- unique(pbmc_obj$cell_type)
  print(paste("🔍 合并前 PBMC 包含类型:", paste(old_types, collapse = ", ")))
  
  # 使用 ifelse 或 recode 进行合并
  # 注意：ScType 数据库中的名称通常为 "Classical Monocytes" 和 "Non-classical monocytes"
  pbmc_obj$cell_type <- ifelse(
    pbmc_obj$cell_type %in% c("Classical Monocytes", "Non-classical monocytes"), 
    "Monocytes", 
    pbmc_obj$cell_type
  )
  
  # 重新转换为 factor 以便后续绘图颜色锁定
  pbmc_obj$cell_type <- factor(pbmc_obj$cell_type)
  
  # 放回列表
  sc_by_tissue[["PBMC"]] <- pbmc_obj
  
  print(paste("✅ PBMC 单核细胞合并完成。当前类型:", paste(unique(pbmc_obj$cell_type), collapse = ", ")))
}

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
# ==============================================================================
# 8. (混合分析版) 深度分析：Fkbp5+ 单核细胞溯源 (Cold_4C vs RT_25C)
# ==============================================================================
print("🚀 步骤8/8: 正在进行 Fkbp5+ 单核细胞的谱系定位分析 (Cold + RT 混合)...")

# ------------------------------------------------------------------------------
# 8.1 提取细胞 (保留 Cold_4C 和 RT_25C)
# ------------------------------------------------------------------------------
myeloid_list <- list()

for (tissue in names(sc_by_tissue)) {
  obj <- sc_by_tissue[[tissue]]
  
  # 1. 查找髓系细胞类型
  target_cells <- grep("Monocytes|Macrophage", obj$cell_type, value = TRUE, ignore.case = TRUE)
  
  if (length(target_cells) > 0) {
    # 2. 提取髓系细胞
    # 【关键修改】条件改为: 属于髓系 且 (属于 Cold 组 或 RT 组)
    sub_obj <- subset(obj, subset = cell_type %in% target_cells & Group %in% c("Cold_4C", "RT_25C"))
    
    if (ncol(sub_obj) > 0) {
      print(paste("  -> 从", tissue, "提取 Cold+RT 髓系细胞:", ncol(sub_obj), "个"))
      sub_obj$Original_Tissue <- tissue
      myeloid_list[[tissue]] <- sub_obj
    }
  }
}

if (length(myeloid_list) == 0) stop("❌ 未找到符合条件的髓系细胞！")

# 合并对象
myeloid_combined <- merge(myeloid_list[[1]], y = myeloid_list[2:length(myeloid_list)])
DefaultAssay(myeloid_combined) <- "RNA"

# ------------------------------------------------------------------------------
# 8.2 数据层合并与 Fkbp5 提取 (Seurat V5 必需)
# ------------------------------------------------------------------------------
print("  -> [V5 修复] 正在合并不同组织的矩阵层 (JoinLayers)...")
myeloid_combined <- JoinLayers(myeloid_combined) 

# 提取 Fkbp5 表达量
if (!"Fkbp5" %in% rownames(myeloid_combined)) {
  real_name <- grep("Fkbp5", rownames(myeloid_combined), ignore.case = TRUE, value = TRUE)
  if(length(real_name) > 0) target_gene <- real_name[1] else stop("❌ 未找到 Fkbp5 基因")
} else {
  target_gene <- "Fkbp5"
}

fkbp5_counts <- GetAssayData(myeloid_combined, assay = "RNA", layer = "counts")[target_gene, ]

# ------------------------------------------------------------------------------
# 8.3 定义详细亚群 (6种可能的组合)
# ------------------------------------------------------------------------------
# 逻辑：先判断是单核还是巨噬 -> 再判断是 Cold 还是 RT -> 再判断(如果是单核)是否 Fkbp5+

cell_types <- myeloid_combined$cell_type
groups <- myeloid_combined$Group
new_labels <- vector("character", length = ncol(myeloid_combined))

for (i in 1:ncol(myeloid_combined)) {
  ctype <- cell_types[i]
  grp   <- groups[i]
  expr  <- fkbp5_counts[i]
  
  # 简化组名后缀
  suffix <- ifelse(grp == "Cold_4C", "(Cold)", "(RT)")
  
  if (grepl("Macrophage", ctype, ignore.case = TRUE)) {
    # 巨噬细胞
    new_labels[i] <- paste("Macrophages", suffix)
  } else {
    # 单核细胞
    if (expr > 0) {
      new_labels[i] <- paste("Fkbp5+ Mono", suffix)
    } else {
      new_labels[i] <- paste("Monocytes", suffix)
    }
  }
}

# 设置因子水平，保证图例顺序清晰
# 顺序：RT单核 -> Cold单核 -> RT Fkbp5+ -> Cold Fkbp5+ (红) -> RT 巨噬 -> Cold 巨噬
my_levels <- c("Monocytes (RT)", "Monocytes (Cold)", 
               "Fkbp5+ Mono (RT)", "Fkbp5+ Mono (Cold)", 
               "Macrophages (RT)", "Macrophages (Cold)")

# 只保留实际存在的 level
existing_levels <- intersect(my_levels, unique(new_labels))
myeloid_combined$Myeloid_Subtype <- factor(new_labels, levels = existing_levels)

print("✅ 分组定义完成。各组数量：")
print(table(myeloid_combined$Myeloid_Subtype))

# ------------------------------------------------------------------------------
# 8.4 重新处理 (标准化 -> 降维)
# ------------------------------------------------------------------------------
print("  -> 正在进行标准化和降维...")
myeloid_combined <- NormalizeData(myeloid_combined)
myeloid_combined <- FindVariableFeatures(myeloid_combined, nfeatures = 2000)
myeloid_combined <- ScaleData(myeloid_combined)
myeloid_combined <- RunPCA(myeloid_combined, verbose = FALSE)

# 此时可以稍微增加 dims，因为包含了更多异质性
myeloid_combined <- RunUMAP(myeloid_combined, dims = 1:25) 

# ------------------------------------------------------------------------------
# 8.5 (优化版) 绘图：三张图独立输出
# ------------------------------------------------------------------------------
print("  -> 正在绘制并保存独立图片...")

# 定义颜色方案 (保持一致性)
color_map <- c(
  "Monocytes (RT)"      = "#D3D3D3",  # LightGrey
  "Monocytes (Cold)"    = "#808080",  # DarkGrey
  "Fkbp5+ Mono (RT)"    = "#FBB4AE",  # Pink
  "Fkbp5+ Mono (Cold)"  = "#E41A1C",  # Red
  "Macrophages (RT)"    = "#A6CEE3",  # LightBlue
  "Macrophages (Cold)"  = "#1F78B4"   # DarkBlue
)

# --- 图 1: 总体分布图 (Lineage Tracing) ---
p_trace <- DimPlot(myeloid_combined, reduction = "umap", group.by = "Myeloid_Subtype", pt.size = 1.2) +
  scale_color_manual(values = color_map) +
  ggtitle("Myeloid Lineage Tracing (Cold vs RT)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

ggsave(filename = "1_Lineage_Tracing_Map.png", plot = p_trace, width = 10, height = 8, path = data_dir)
print("  ✅ [1/3] 总体分布图已保存: 1_Lineage_Tracing_Map.png")

# --- 图 2: Cold_4C 下的 Fkbp5 表达 ---
# 提取 Cold 细胞
obj_cold <- subset(myeloid_combined, subset = Group == "Cold_4C")

p_cold_gene <- FeaturePlot(obj_cold, features = target_gene, order = TRUE, pt.size = 1.2) + 
  scale_color_viridis_c(option = "plasma") + # 使用对比度更高的 plasma 色盘
  ggtitle(paste0(target_gene, " Expression (Cold_4C)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(filename = "2_Fkbp5_Expression_Cold.png", plot = p_cold_gene, width = 8, height = 7, path = data_dir)
print("  ✅ [2/3] 寒冷组基因图已保存: 2_Fkbp5_Expression_Cold.png")

# --- 图 3: RT_25C 下的 Fkbp5 表达 ---
# 提取 RT 细胞
obj_rt <- subset(myeloid_combined, subset = Group == "RT_25C")

p_rt_gene <- FeaturePlot(obj_rt, features = target_gene, order = TRUE, pt.size = 1.2) + 
  scale_color_viridis_c(option = "plasma") + 
  ggtitle(paste0(target_gene, " Expression (RT_25C)")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(filename = "3_Fkbp5_Expression_RT.png", plot = p_rt_gene, width = 8, height = 7, path = data_dir)
print("  ✅ [3/3] 常温组基因图已保存: 3_Fkbp5_Expression_RT.png")

print("🎉 所有图片绘制完成！")
