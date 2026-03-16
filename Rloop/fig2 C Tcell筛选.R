# 分别运行 ProjecTILs 核心对齐与注释功能
CD4_obj <- Run.ProjecTILs(CD4_obj, ref = ref_CD4)
CD8_obj <- Run.ProjecTILs(CD8_obj, ref = ref_CD8)

# 提取并合并注释结果 (functional.cluster)
cd4_meta <- CD4_obj@meta.data[, "functional.cluster", drop = FALSE]
cd8_meta <- CD8_obj@meta.data[, "functional.cluster", drop = FALSE]
combined_meta <- rbind(cd4_meta, cd8_meta)

# 将合并后的 functional.cluster 列重新添加回总的 Tcell 对象
Tcell <- AddMetaData(Tcell, metadata = combined_meta)

# 对于少部分既不表达 CD4 也不表达 CD8A 的双阴性细胞，标记为普通 T 细胞
Tcell$functional.cluster[is.na(Tcell$functional.cluster)] <- "Unclassified_T"
# --- 修改结束 ---

# 6. 将 ProjecTILs 的标准命名转换为你论文图纸里的标签风格
# ProjecTILs 会生成一列叫 functional.cluster 的结果
# 我们将其重命名并映射，存入一个新的列 `Tcell_subtype`
Tcell@meta.data <- Tcell@meta.data %>%
  mutate(Tcell_subtype = case_when(
    functional.cluster == "Cd8_Tex" ~ "CD8 Exhaust",
    functional.cluster == "Cd8_Teff" ~ "CD8 Cytotoxic",
    functional.cluster == "Cd8_NaiveLike" ~ "CD8 Naive",
    functional.cluster == "Cd8_Memory" ~ "CD8 Memory", # 图中未展示，但图谱有
    functional.cluster == "Cd4_Treg" ~ "CD4 Treg",
    functional.cluster == "Cd4_NaiveLike" ~ "CD4 Naive",
    functional.cluster == "Th1" ~ "CD4 Th1",
    functional.cluster == "Tfh" ~ "CD4 Tfh",
    TRUE ~ as.character(functional.cluster) # 保留其他未匹配的原始名称
  ))

# 为了适配你后续的 C 部分画图代码，我们将细分结果也同步给 singleR
# 这样你跑后续堆叠柱状图时，就能直接出来各个亚群的比例了
Tcell@meta.data$singleR <- Tcell@meta.data$Tcell_subtype

# 7. 打印结果检查
print("T 细胞亚群分类结果统计：")
print(table(Tcell@meta.data$singleR))

# 8. 覆盖保存为 Tcell.Rdata
save(Tcell, file = "Tcell.Rdata")
