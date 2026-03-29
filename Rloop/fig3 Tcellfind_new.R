setwd("/mnt/disk1/qiuzerui/downloads/Rloop")
load("scRNA.Rdata")
# 1. 提取 SingleR 注释为 T_cells 的细胞
t_cells_subset <- subset(scRNA, subset = singleR == "T_cells")
DefaultAssay(t_cells_subset) <- 'RNA'
# 2. 重新运行标准流程 (非常重要：为了重新寻找 T 细胞内部的高变基因)
# 如果你之前用的是标准流程：
t_cells_subset <- NormalizeData(t_cells_subset)
t_cells_subset <- FindVariableFeatures(t_cells_subset, nfeatures = 2000)
t_cells_subset <- ScaleData(t_cells_subset)
t_cells_subset <- RunPCA(t_cells_subset)

# 3. 重新降维聚类 (建议 dims 可以设小一点，比如 10-15)
t_cells_subset <- RunUMAP(t_cells_subset, dims = 1:15)
t_cells_subset <- FindNeighbors(t_cells_subset, dims = 1:15)
t_cells_subset <- FindClusters(t_cells_subset, resolution = 0.5) # 分辨率可调
DimPlot(scRNA, group.by="seurat_clusters", label=FALSE, label.size=5, reduction='umap')
# ==============================================================================
# 4. 使用 AddModuleScore 进行 T 细胞状态评分 (保持不变)
# ==============================================================================
t_features_list <- list(
  Naive = c("CCR7", "LEF1", "SELL", "TCF7", "NOSIP"),
  Cytotoxic = c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "FGFBP2"),
  Exhausted = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "LAYN", "TOX"),
  Treg = c("FOXP3", "IL2RA", "IKZF2", "BATF", "LAIR2"),
  Th17 = c("RORC", "IL17A", "IL17F", "CCR6", "KLRB1")
)

t_cells_subset <- AddModuleScore(
  object = t_cells_subset,
  features = t_features_list,
  ctrl = 100,
  name = "T_State_"
)

score_names <- names(t_features_list)
for(i in 1:length(score_names)){
  colnames(t_cells_subset@meta.data)[which(colnames(t_cells_subset@meta.data) == paste0("T_State_", i))] <- paste0(score_names[i], "_Score")
}

# ==============================================================================
# 5. 亚群辅助鉴定与可视化 (增加谱系检查)
# ==============================================================================
FeaturePlot(t_cells_subset, features = c("CD4", "CD8A"), ncol = 2)

# ==============================================================================
# 6. 核心修改：先区分 CD4/CD8 谱系，再合并亚群状态
# ==============================================================================

# 1. 获取 CD4 和 CD8A 的标准化表达量
# 使用 Seurat v5 的 LayerData 提取数据
cd4_expr <- LayerData(t_cells_subset, assay = "RNA", layer = "data")["CD4", ]
cd8_expr <- LayerData(t_cells_subset, assay = "RNA", layer = "data")["CD8A", ]

# 2. 判定谱系 (Lineage): 比较 CD4 和 CD8A 的表达高低
t_cells_subset$lineage <- ifelse(cd4_expr >= cd8_expr, "CD4", "CD8")

# 3. 判定状态 (State): 基于 Score 最大值 (原有逻辑)
scores_matrix <- t_cells_subset@meta.data[, paste0(score_names, "_Score")]
t_cells_subset$predicted_state <- score_names[max.col(scores_matrix)]

# 4. 复合命名：将谱系与状态拼接，并根据生物学常规修正 (如 Treg/Th17 默认归为 CD4)
t_cells_subset$celltype_fine <- paste(t_cells_subset$lineage, t_cells_subset$predicted_state)

# 修正：根据你图片中的分类，Treg 和 Th17 统一前缀为 CD4
t_cells_subset$celltype_fine <- gsub("CD8 Treg", "CD4 Treg", t_cells_subset$celltype_fine)
t_cells_subset$celltype_fine <- gsub("CD8 Th17", "CD4 Th17", t_cells_subset$celltype_fine)

# 设置为默认 Ident 方便绘图
Idents(t_cells_subset) <- "celltype_fine"

# 检查最终分类结果
DimPlot(t_cells_subset, reduction = "umap", label = TRUE)

# 最后保存结果
save(t_cells_subset, file = "T_cells_Subsets_with_Scores.Rdata")
# ==============================================================================
# 7. 导出细胞注释结果
# ==============================================================================

# 提取 Barcode 和对应的细分亚群标签
# rownames(t_cells_subset@meta.data) 获取的是细胞唯一的条形码
output_df <- data.frame(
  CellName = rownames(t_cells_subset@meta.data),
  celltype_fine = t_cells_subset$celltype_fine
)

# 导出为 Tab 分隔的文本文件
write.table(output_df, 
            file = "T_cell_fine_annotation.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# 在控制台打印前几行确认一下
print("导出成功！前 5 行预览：")
print(head(output_df))
