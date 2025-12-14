# ==============================================================================
# 单细胞转录组全流程分析代码 (Singleron/10x 格式) - 最终修复版
# 实验设计: 
#   - 组织: Aorta (A), PBMC (B), BoneMarrow (M)
#   - 分组: 1=25°C (RT), 2=4°C (Cold), 3=30°C (TN)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx) # 用于读取 ScType 数据库
library(Matrix)   # 用于处理稀疏矩阵
library(HGNChelper)

# 加载 ScType 核心函数 (需保持联网)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# ------------------------------------------------------------------------------
# 1. 设置路径与元数据映射
# ------------------------------------------------------------------------------
# !!! 请务必修改为您电脑上的实际路径 (注意使用 / 而不是 \) !!!
data_dir <- "D:/qiuzerui/单细胞寒冷处理小鼠" 

# 样本ID列表
sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")

# 组织映射: A=主动脉, B=外周血, M=骨髓
tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")

# 分组映射: 1=25度(室温), 2=4度(寒冷), 3=30度(热中性)
group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")

# ------------------------------------------------------------------------------
# 2. 读取数据、修复格式并添加 Metadata
# ------------------------------------------------------------------------------
sc_list <- list()
full_folder_names <- list.files(data_dir, pattern = "_matrix")

print("🚀 步骤1/5: 开始读取并修复 Singleron 数据...")

for (sample in sample_ids) {
  # --- A. 匹配文件夹 ---
  folder <- full_folder_names[grep(paste0("^", sample, "_"), full_folder_names)]
  
  if(length(folder) == 0) {
    message(paste("⚠️ 跳过：未找到样本", sample)); next
  }
  
  # --- B. 读取原始数据 ---
  print(paste("   正在处理:", sample))
  
  # 使用 tryCatch 捕获读取错误
  tryCatch({
    counts <- Read10X(data.dir = file.path(data_dir, folder[1]))
    
    # 兼容性修复：如果 Read10X 返回的是列表，提取 Gene Expression 矩阵
    if (is.list(counts) && !is(counts, "dgCMatrix")) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts$`Gene Expression` 
      } else {
        counts <- counts[[1]]
      }
    }
    
    # --- C. 关键修复：解决 "No cell overlap" 报错 (1/2) ---
    # 将 Singleron 细胞名中的下划线 "_" 替换为 Seurat 喜欢的减号 "-"
    colnames(counts) <- gsub("_", "-", colnames(counts))
    
    # --- D. 创建对象 ---
    # min.cells=3: 去除极少出现的基因
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    # 检查是否因为过滤导致细胞数为0
    if (ncol(sc_obj) == 0) {
      message(paste("   ⚠️ 警告：样本", sample, "创建后无细胞，跳过。")); next
    }
    
    # --- E. 添加 Metadata (关键修复点 2/2) ---
    prefix <- substr(sample, 1, 1) # A, B, M
    suffix <- substr(sample, 2, 2) # 1, 2, 3
    
    # !!! 这里的 as.character 是修复 "No cell overlap" 的核心 !!!
    # 它去掉了 tissue_map 向量的名字 "A"，只保留值 "Aorta"
    sc_obj$Tissue <- as.character(tissue_map[prefix])
    sc_obj$Group  <- as.character(group_map[suffix])
    
    # 标记 Condition (方便后续差异分析)
    if(suffix == "2") {
      sc_obj$Condition <- "Cold"
    } else {
      sc_obj$Condition <- "NonCold"
    }
    
    # 计算线粒体比例
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
    
    sc_list[[sample]] <- sc_obj
    print(paste("   ✅ 成功入库: ", sample, "| 细胞数:", ncol(sc_obj), "| 组织:", sc_obj$Tissue[1]))
    
  }, error = function(e) {
    message(paste("   ❌ 处理样本", sample, "时出错:", e$message))
  })
}

# 合并所有样本
if (length(sc_list) > 0) {
  print("正在合并所有样本...")
  sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = sample_ids)
} else {
  stop("❌ 未读取到任何数据，请检查路径设置！")
}

# ------------------------------------------------------------------------------
# 3. 统一质控 (QC)
# ------------------------------------------------------------------------------
print("🚀 步骤2/5: 正在进行质控过滤...")
# 过滤标准：
# 1. nFeature_RNA > 200: 去除死细胞/碎片
# 2. nFeature_RNA < 6000: 去除双细胞
# 3. percent.mt < 15: 去除线粒体过高的濒死细胞 (Aorta样本如果质量差可适当放宽到20)
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# ------------------------------------------------------------------------------
# 4. 按组织拆分并独立分析 (聚类 + 降维)
# ------------------------------------------------------------------------------
print("🚀 步骤3/5: 按组织拆分并独立聚类...")
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

# 定义标准处理流程函数
run_standard_pipeline <- function(obj, tissue_name) {
  print(paste("   正在处理组织:", tissue_name, "..."))
  
  # 简单检查细胞数
  if(ncol(obj) < 50) { return(NULL) }
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5) # Resolution 越大，分群越细
  obj <- RunUMAP(obj, dims = 1:20)
  return(obj)
}

# 批量运行
sc_by_tissue <- lapply(names(sc_by_tissue), function(x) {
  run_standard_pipeline(sc_by_tissue[[x]], x)
})
names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow") 
# 移除空对象
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]

# ------------------------------------------------------------------------------
# 5. ScType 自动细胞注释 (本地文件版)
# ------------------------------------------------------------------------------
print("🚀 步骤4/5: 运行 ScType 细胞注释...")

# !!! 修改这里：不再使用 https 链接，而是使用本地文件路径 !!!
# 确保你已经把下载的 ScTypeDB_full.xlsx 放到了 data_dir 里
db_file_path <- file.path(data_dir, "ScTypeDB_full.xlsx") 

# 检查文件是否存在
if (!file.exists(db_file_path)) {
  stop(paste("❌ 未找到数据库文件！请手动下载 ScTypeDB_full.xlsx 并放入", data_dir))
}

# 准备数据库 (针对免疫系统)
# gene_sets_prepare 函数可以直接读取本地 xlsx 文件
gs_list_immune <- gene_sets_prepare(db_file_path, "Immune system") 

# 注释函数 (保持不变)
run_annotation <- function(obj, gs_list) {
  # 计算打分
  es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(obj, layer="scale.data")), scaled = TRUE, 
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # 将打分映射到 Cluster
  cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # 写入 Metadata
  obj@meta.data$cell_type <- ""
  for(j in unique(sctype_scores$cluster)){
    cl_type <- sctype_scores[sctype_scores$cluster==j, "type"]
    obj@meta.data$cell_type[obj@meta.data$seurat_clusters == j] <- as.character(cl_type)
  }
  return(obj)
}

# 对三个组织分别注释
sc_by_tissue <- lapply(sc_by_tissue, function(x) run_annotation(x, gs_list_immune))


# ==============================================================================
# 批量绘制并保存所有组织 (Aorta, PBMC, BoneMarrow) 的注释结果
# ==============================================================================


# 检查结果列表里都有哪些组织
print(paste("当前包含的组织有:", paste(names(sc_by_tissue), collapse = ", ")))

# 循环遍历每个组织进行绘图
for (tissue_name in names(sc_by_tissue)) {
  
  print(paste("正在绘制:", tissue_name, "..."))
  
  # 1. 获取当前组织的对象
  obj <- sc_by_tissue[[tissue_name]]
  
  # 2. 生成 UMAP 图 (按 cell_type 着色)
  p <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "- ScType Annotation")) +
    theme(plot.title = element_text(hjust = 0.5)) # 标题居中
  
  # 3. 定义文件名
  #filename_pdf <- paste0(tissue_name, "_Annotation.pdf")
  filename_png <- paste0(tissue_name, "_Annotation.png")
  
  # 4. 保存图片 (避免窗口报错，直接存文件)
  # 这里的 data_dir 是你之前设置的文件夹路径
  tryCatch({
    #ggsave(filename = filename_pdf, plot = p, width = 10, height = 8, path = data_dir)
    ggsave(filename = filename_png, plot = p, width = 10, height = 8, path = data_dir)
    print(paste("   ✅ 图片已保存:", filename_png))
  }, error = function(e) {
    print(paste("   ❌ 保存失败:", e$message))
  })
}

print("🎉 所有组织的图片绘制完成！请去文件夹查看。")

# ------------------------------------------------------------------------------
# 6. 差异表达分析 (Differential Expression)
# ------------------------------------------------------------------------------
print("🚀 步骤5/5: 差异表达分析 (Cold vs RT / Cold vs TN)...")

# 定义差异分析函数
analyze_differences <- function(tissue_name, target_cell_type) {
  
  if (!tissue_name %in% names(sc_by_tissue)) return(NULL)
  obj <- sc_by_tissue[[tissue_name]]
  
  # 检查是否存在该细胞
  if (!target_cell_type %in% obj$cell_type) {
    message(paste("   ⚠️ 跳过: 在", tissue_name, "中未找到", target_cell_type)); return(NULL)
  }
  
  print(paste("   正在分析:", tissue_name, "-", target_cell_type))
  
  # 1. 提取特定细胞
  Idents(obj) <- "cell_type"
  sub_obj <- subset(obj, idents = target_cell_type)
  
  # 2. 切换到分组比较
  Idents(sub_obj) <- "Group"
  
  # --- 比较 A: 寒冷(2) vs 室温(1) ---
  # logFC > 0 表示在寒冷组上调
  try({
    de_2_vs_1 <- FindMarkers(sub_obj, ident.1 = "Cold_4C", ident.2 = "RT_25C", min.pct = 0.25)
    write.csv(de_2_vs_1, paste0(tissue_name, "_", gsub(" ", "_", target_cell_type), "_Cold_vs_RT.csv"))
  }, silent=TRUE)
  
  # --- 比较 B: 寒冷(2) vs 热中性(3) ---
  # logFC > 0 表示寒冷组比热中性组高
  try({
    de_2_vs_3 <- FindMarkers(sub_obj, ident.1 = "Cold_4C", ident.2 = "TN_30C", min.pct = 0.25)
    write.csv(de_2_vs_3, paste0(tissue_name, "_", gsub(" ", "_", target_cell_type), "_Cold_vs_TN.csv"))
  }, silent=TRUE)
}

# === 运行示例 (您可以根据 ScType 的结果修改下面的细胞名字) ===

# ==============================================================================
# 针对你提供的细胞类型进行差异分析
# ==============================================================================
# ------------------------------------------------------------------------------
# 6. 精准靶向差异分析 (B组->T细胞, A组->巨噬, M组->HSC)
#请修改参数
# ------------------------------------------------------------------------------
print("🚀 步骤5/6: 核心修复 - 合并数据层 (JoinLayers)...")
for (tissue in names(sc_by_tissue)) {
  sc_by_tissue[[tissue]] <- JoinLayers(sc_by_tissue[[tissue]])
}

print("🚀 步骤6/6: 开始一对一精准靶向分析...")

# 定义差异分析核心函数
analyze_differences <- function(tissue_name, target_cell_type) {
  obj <- sc_by_tissue[[tissue_name]]
  
  # 提取细胞
  Idents(obj) <- "cell_type"
  sub_obj <- subset(obj, idents = target_cell_type)
  Idents(sub_obj) <- "Group"
  groups_here <- unique(sub_obj$Group)
  
  print(paste("   📊 分析:", tissue_name, "->", target_cell_type))
  
  # Cold vs RT
  if ("Cold_4C" %in% groups_here && "RT_25C" %in% groups_here) {
    try({
      markers <- FindMarkers(sub_obj, ident.1 = "Cold_4C", ident.2 = "RT_25C", min.pct = 0.1, logfc.threshold = 0.25)
      markers$gene <- rownames(markers)
      file_name <- paste0(tissue_name, "_", gsub("[ /+]", "_", target_cell_type), "_Cold_vs_RT.csv")
      write.csv(markers, file.path(data_dir, file_name), row.names = FALSE)
      print(paste("      ✅ 保存:", file_name))
    }, silent = FALSE)
  }
}

# 定义您的靶向任务
tasks <- list(
  list(tissue = "PBMC",       keyword = "T cell"),          # B组找 T细胞
  list(tissue = "Aorta",      keyword = "Macrophage"),      # A组找 巨噬细胞
  list(tissue = "BoneMarrow", keyword = "Progenitor|HSC")   # M组找 HSC
)

# 执行任务
for (task in tasks) {
  target_tissue <- task$tissue
  target_key    <- task$keyword
  
  if (!target_tissue %in% names(sc_by_tissue)) next
  
  # 模糊搜索匹配的细胞
  all_cells <- unique(sc_by_tissue[[target_tissue]]$cell_type)
  matched_cells <- grep(target_key, all_cells, value = TRUE, ignore.case = TRUE)
  
  if (length(matched_cells) > 0) {
    print(paste("📂 [", target_tissue, "] 匹配到:", paste(matched_cells, collapse=", ")))
    for (cell in matched_cells) {
      analyze_differences(target_tissue, cell)
    }
  } else {
    print(paste("⚠️ [", target_tissue, "] 未找到关键词:", target_key))
  }
}

print("🎉 全流程运行完毕！结果已保存在您的文件夹中。")
#7.拟时序分析
# ==============================================================================
# PBMC T细胞拟时序分析 (独立版 - 不依赖 SeuratWrappers)
# 解决 GitHub 安装报错的最佳方案
# ==============================================================================

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)

# ------------------------------------------------------------------------------
# 1. 准备数据：提取 PBMC 中的 T 细胞
# ------------------------------------------------------------------------------
print("⏳ 正在从 PBMC 中提取 T 细胞...")

# 检查对象是否存在
if (!"PBMC" %in% names(sc_by_tissue)) stop("❌ 错误：找不到 PBMC 对象，请先运行前面的步骤。")

pbmc <- sc_by_tissue[["PBMC"]]
Idents(pbmc) <- "cell_type"

# 模糊匹配 T 细胞 (自动找 Naive, Memory, CD8 等)
t_cell_names <- grep("T cell", unique(pbmc$cell_type), value = TRUE, ignore.case = TRUE)

if (length(t_cell_names) == 0) stop("❌ 未在 PBMC 中找到 T 细胞！")

print(paste("   发现 T 细胞亚群:", paste(t_cell_names, collapse = ", ")))
sub_obj <- subset(pbmc, idents = t_cell_names)

# ------------------------------------------------------------------------------
# 2. 手动构建 Monocle3 对象 (关键步骤：替代 SeuratWrappers)
# ------------------------------------------------------------------------------
print("⏳ 正在手动构建 Monocle3 对象...")

# A. 提取表达矩阵 (Counts)
# Seurat V5 需要用 layer="counts"
# 如果报错，说明是旧版 Seurat，会自动切换到 "counts" slot
data_matrix <- tryCatch({
  GetAssayData(sub_obj, layer = "counts")
}, error = function(e) {
  GetAssayData(sub_obj, slot = "counts")
})

# B. 提取细胞元数据 (Metadata)
cell_metadata <- sub_obj@meta.data

# C. 提取基因元数据 (Gene Metadata)
gene_metadata <- data.frame(gene_short_name = rownames(data_matrix))
rownames(gene_metadata) <- rownames(data_matrix)

# D. 创建 CDS 对象
cds <- new_cell_data_set(data_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# ------------------------------------------------------------------------------
# 3. 轨迹构建流程
# ------------------------------------------------------------------------------
print("⏳ 正在预处理数据...")

# A. 预处理 (PCA)
cds <- preprocess_cds(cds, num_dim = 30)

# B. 关键点：强制使用 Seurat 的 UMAP 坐标
# 这样 Monocle 画出来的图和 Seurat 的图长得一模一样
print("   同步 Seurat UMAP 坐标...")
seurat_umap <- sub_obj@reductions$umap@cell.embeddings
# 确保细胞顺序一致
seurat_umap <- seurat_umap[colnames(cds), ]
# 覆盖 Monocle 的坐标
cds@int_colData@listData$reducedDims$UMAP <- seurat_umap

# C. 聚类 (Monocle 必需步骤，用于划分轨迹)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# D. 学习轨迹 (Learn Graph)
print("⏳ 正在构建轨迹图 (Learn Graph)...")
cds <- learn_graph(cds, use_partition = TRUE)

# ------------------------------------------------------------------------------
# 4. 自动寻找起点 (Naive + RT) 并定根
# ------------------------------------------------------------------------------
print("⏳ 正在计算伪时间 (Pseudotime)...")

get_t_cell_root <- function(cds) {
  meta <- colData(cds)
  # 优先找: 名字带 "Naive" 且分组是 "RT_25C" 的细胞
  candidate_cells <- rownames(meta)[grepl("Naive", meta$cell_type, ignore.case = TRUE) & 
                                      meta$Group == "RT_25C"]
  
  # 如果找不到 (比如 ScType 没注释出 Naive)，退而求其次：找所有 RT_25C 细胞
  if (length(candidate_cells) < 5) {
    print("   (提示: Naive 细胞不足，放宽条件至所有 RT 细胞)")
    candidate_cells <- rownames(meta)[meta$Group == "RT_25C"]
  }
  
  # 找到这些细胞在 UMAP 上最密集的那个节点
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[candidate_cells,]))))]
  return(root_node)
}

# 尝试自动定根
tryCatch({
  root_node <- get_t_cell_root(cds)
  cds <- order_cells(cds, root_pr_nodes = root_node)
  print("✅ 伪时间计算完成！")
}, error = function(e) {
  print("⚠️ 自动定根失败 (可能是细胞太少)，跳过定根步骤。")
})

# ------------------------------------------------------------------------------
# 5. 绘图与保存
# ------------------------------------------------------------------------------
print("🎨 正在绘图...")

# 图 A: 伪时间 (颜色越亮 = 时间越晚)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                 graph_label_size=1.5, cell_size = 0.8) + ggtitle("Pseudotime (Time)")

# 图 B: 分组 (看红色的 Cold 组是不是在末端)
p2 <- plot_cells(cds, color_cells_by = "Group", 
                 label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                 graph_label_size=1.5, cell_size = 0.8) + ggtitle("Group Distribution")

# 图 C: 细胞类型 (验证起点是不是 Naive)
p3 <- plot_cells(cds, color_cells_by = "cell_type", 
                 label_cell_groups=TRUE, label_leaves=FALSE, label_branch_points=FALSE,
                 graph_label_size=3, cell_size = 0.8) + ggtitle("Cell Types")

# 拼图并保存
combined_plot <- (p1 | p2) / p3
ggsave(file.path(data_dir, "PBMC_Tcells_Trajectory_Final.png"), combined_plot, width = 12, height = 12)

print(paste("🎉 分析完成！图片已保存至:", file.path(data_dir, "PBMC_Tcells_Trajectory_Final.png")))
