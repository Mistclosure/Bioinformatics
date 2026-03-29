#D----
library(GSVA)
library(Seurat) # 确保加载 Seurat
library(GSEABase)
load("scRNA.Rdata")

# Seurat v5：合并层级，确保 SCT 或 RNA 层的数据完整性
scRNA <- JoinLayers(scRNA)

meta.data = scRNA@meta.data
meta.data = meta.data[meta.data$singleR == "NK_cell",]
B = as.data.frame(table(meta.data$orig.ident))

# 提取 NK 细胞子集
scRNA = scRNA[, rownames(meta.data)]

geneSets <- getGmt('GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP.v2026.1.Hs.gmt')

# --- Seurat v5 更新：使用 LayerData 提取 SCT 层的 data ---
# 确保 SCT assay 是激活状态
DefaultAssay(scRNA) <- "SCT"
exp = as.data.frame(LayerData(scRNA, assay = "SCT", layer = "data"))

dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

# 运行 GSVA (ssgsea)
# 1. 构造 ssgsea 的参数对象
# 注意：新版本中 exprData 对应你的矩阵，geneSets 对应你的 GMT
# ssgsea 不需要 mx.diff 参数，那个是给标准 GSVA 方法用的
ssgsea_param <- ssgseaParam(exprData = mat, 
                            geneSets = geneSets,
                            minSize = 1,       # 基因集的最小长度
                            maxSize = Inf,     # 基因集的最大长度
                            normalize = TRUE)  # 是否进行归一化

# 2. 运行 gsva，直接传入参数对象即可
GSVA_hall <- gsva(ssgsea_param)

# --- 修复变量名：统一使用 GSVA_hall ---
hall = t(GSVA_hall) 
hall = as.data.frame(hall)
# 这里的 meta.data[1] 提取的是 orig.ident 所在的列
hall = merge(meta.data[1], hall, by=0)
hall_sample <- hall %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
Type = read.table("Rloopscore.txt", header=T, sep="\t", check.names=F, row.names = 1)

# 合并 Rloop 分组信息
rt = merge(Type, hall_sample, by.x = 1, by.y = 1) # 修复 merge 逻辑确保 ID 匹配
rt3 = rt[, c(3, 4)]
colnames(rt3) = c("Type", "NK & T-cell exhaustion")

group=levels(factor(rt3$Type))
rt3$Type=factor(rt3$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

# 绘图
ggplot(rt3, aes(x=Type, y=`NK & T-cell exhaustion`, fill=Type)) +
  geom_violin(trim=FALSE, color="white") +
  geom_boxplot(width=0.2, position=position_dodge(0.9)) +
  geom_point(aes(x=Type, y=`NK & T-cell exhaustion`), pch=19, position=position_dodge(0.9), size=0.5) +
  scale_fill_manual(values = c("#eaeae8","#7bacb0")) +
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons, method="t.test") +
  theme(
    legend.position = "none",
    panel.border = element_blank(), 
    axis.line = element_line(colour = "black", size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("") +
  ylab("NK exhaustion Score") + xlab("")
