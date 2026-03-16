
#首先二次注释Tcells
# 1. 加载必要的 R 包
library(Seurat)
library(ProjecTILs)
library(dplyr)

# 2. 加载全量单细胞数据
# 请确保此时你在正确的工作目录下
load("scRNA.Rdata")

# 3. 提取 T 细胞子集
# 根据你之前的代码，T 细胞在 singleR 列中标记为 'T_cells'
# 如果有其他名称（如 'CD8+ T' 等），可以用 %in% c("T_cells", "...") 组合
Tcell <- subset(scRNA, subset = singleR == "T_cells")

# 4. 加载 ProjecTILs 默认的肿瘤浸润 T 细胞参考图谱
ref <- load.reference.map(ref = "Human_STACAS_v1")

# 5. 运行 ProjecTILs 核心对齐与注释功能
# ProjecTILs 会自动将你的人类基因 (大写) 映射到其参考图谱中
Tcell <- Run.ProjecTILs(Tcell, ref = ref)

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
