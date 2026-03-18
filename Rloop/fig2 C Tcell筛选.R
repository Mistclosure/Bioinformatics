# 1. 加载必要的 R 包
library(Seurat)
library(ProjecTILs)
library(dplyr)
setwd('/mnt/disk1/qiuzerui/downloads/Rloop/')
# 2. 加载全量单细胞数据
load("scRNA.Rdata")

# 3. 提取 T 细胞子集
Tcell <- subset(scRNA, subset = singleR == "T_cells")

# --- 修改开始：群划分逻辑与精细注释 ---

# 4. 加载针对人类 CD4 和 CD8 的特定参考图谱 (读取本地文件)
ref_CD4 <- readRDS("CD4T_human_ref_v1.rds")
ref_CD8 <- readRDS("CD8T_human_ref_v1.rds")

# 5. 基于【群】划分 CD4 和 CD8 阵营 (解决 Dropout 导致的 Unclassified 问题)
DefaultAssay(Tcell) <- "RNA"
Tcell <- JoinLayers(Tcell)
Tcell <- NormalizeData(Tcell)

# A. 计算每个 Cluster 的标记基因平均值 (同时考虑 CD8A 和 CD8B)
cluster_averages <- AverageExpression(Tcell, features = c("CD4", "CD8A", "CD8B"), group.by = "seurat_clusters")$RNA

# B. 判定每个 Cluster 的归属
cluster_lineage <- apply(cluster_averages, 2, function(x) {
  cd4_val <- x["CD4"]
  cd8_val <- max(x["CD8A"], x["CD8B"]) 
  if (is.na(cd4_val) || is.na(cd8_val)) return("Unknown") # 防止基因缺失
  if (cd4_val > cd8_val) return("CD4") else return("CD8")
})

# --- 修复后的 C 部分 ---

# 【核心修复】：去掉名字里的 'g'
names(cluster_lineage) <- gsub("^g", "", names(cluster_lineage))

# 【关键改动】：使用 unname() 剥离名字，防止 Seurat 误以为 Cluster ID 是 Cell ID
Tcell$T_lineage <- unname(cluster_lineage[as.character(Tcell$seurat_clusters)])

# 验证一下：现在应该能打印出正常的数字了
print("阵营划分统计：")
print(table(Tcell$T_lineage))

# D. 切分并投影 (这里逻辑不变)
CD4_obj <- subset(Tcell, subset = T_lineage == "CD4")
CD8_obj <- subset(Tcell, subset = T_lineage == "CD8")

CD4_obj <- Run.ProjecTILs(CD4_obj, ref = ref_CD4)
CD8_obj <- Run.ProjecTILs(CD8_obj, ref = ref_CD8)

# --- 5. E. 【最终修复版】：直接映射标签，不再使用 AddMetaData ---

# 1. 建立一个“Cell ID -> 标签”的字典向量
labels_vec <- c(
  setNames(as.character(CD4_obj$functional.cluster), rownames(CD4_obj@meta.data)),
  setNames(as.character(CD8_obj$functional.cluster), rownames(CD8_obj@meta.data))
)

# 2. 直接根据 Tcell 的行名（Cell ID）从字典里拔出对应的标签
# unname() 非常关键：它能防止 Seurat 自作聪明地去比对名字，从而跳过报错
Tcell$functional.cluster <- unname(labels_vec[rownames(Tcell@meta.data)])

# 3. 兜底处理：将那些由于双阴性没被分配到 CD4/CD8 阵营的细胞标记
Tcell$functional.cluster[is.na(Tcell$functional.cluster)] <- "Unclassified_T"

# --- 6. 将 ProjecTILs 的命名转换为论文标签风格 (逻辑衔接) ---
# 此时 functional.cluster 已经稳稳地存在于 Tcell@meta.data 里的每一行了
Tcell@meta.data <- Tcell@meta.data %>%
  mutate(Tcell_subtype = case_when(
    # CD8 亚型对应
    functional.cluster == "CD8.TEX" ~ "CD8 Exhaust",
    functional.cluster == "CD8.TPEX" ~ "CD8 Pro-Exhaust",
    functional.cluster == "CD8.TEMRA" ~ "CD8 Cytotoxic",
    functional.cluster == "CD8.EM" ~ "CD8 Memory",
    functional.cluster == "CD8.CM" ~ "CD8 Central Memory",
    functional.cluster == "CD8.NaiveLike" ~ "CD8 Naive",
    # CD4 亚型对应
    functional.cluster == "CD4.CTL_Exh" ~ "CD4 Exhaust",
    functional.cluster == "CD4.Treg" ~ "CD4 Treg",
    functional.cluster == "CD4.Tfh" ~ "CD4 Tfh",
    functional.cluster == "CD4.NaiveLike" ~ "CD4 Naive",
    functional.cluster == "CD4.Th17" ~ "CD4 Th17",
    # 兼容处理
    TRUE ~ as.character(functional.cluster) 
  ))

# 同步给 singleR，方便你后面的饼图脚本直接运行
Tcell$singleR <- Tcell$Tcell_subtype

# 为了适配你后续的 C 部分画图代码
Tcell@meta.data$singleR <- Tcell@meta.data$Tcell_subtype

# 7. 打印结果检查
print("T 细胞亚群分类结果统计：")
print(table(Tcell@meta.data$singleR))

# 8. 覆盖保存为 Tcell.Rdata
save(Tcell, file = "Tcell.Rdata")
