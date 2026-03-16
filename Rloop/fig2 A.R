# 加载必要的 R 包
library(Seurat)
library(ggplot2)
library(GSEABase)
library(GSVA)
library(ggpubr)
library(BiocParallel) # 添加并行运算包

# fig2----
# A----
# 加载预先保存的单细胞数据对象
setwd('/mnt/disk1/qiuzerui/downloads/Rloop/')
load("Malignant.Rdata")

# 读取外部定义的评分文件
score <- read.table("score.txt", sep = "\t", header = T, check.names = F, row.names = 1)

# 以中位数为界，将所有细胞分为 "High" 和 "Low" 两组
score$Type <- ifelse(score$score > median(score$score), "High", "Low")

# 确保评分矩阵的行名与单细胞对象的列名完全匹配
score <- score[colnames(pbmc1), ]

# 将分组信息存入单细胞对象的元数据中
pbmc1$Type <- score$Type

# --- 直接将基因列表写入脚本，替代读取 TAA.txt ---
all_genes <- c(
  # TAA 基因
  "TP53", "KRT19", "GRP", "KRT20", "EGFR", "KRAS", "NAPSA", "KRT7", 
  "AKAP3", "DDX43", "PRAME", "PLAC1", "XAGE2", "CALR3", "HORMAD2", 
  "CAGE1", "PAGE5", "XAGE3", "ACTL8", "SAGE1", "MAGEA6", "MAGEB2", 
  "MAGEA1", "MAGEA10", "MAGEA4", "MAGEA11", "MAGEB6", "MAGEB1", "MAGEA9",
  # MHC 基因
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-DPA1", 
  "HLA-DQB1", "HLA-DQB2", "HLA-DQA2", "HLA-DRB5",
  # Immune suppression 基因
  "TGFBR1", "TGFBR2", "TGFBR3", "IL10", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "CD276"
)

# 使用 Seurat 生成气泡图，直接使用定义的 all_genes
P1 <- DotPlot(pbmc1, features = all_genes, group.by = "Type") + coord_flip()

# 提取气泡图底层数据准备美化
rt <- P1[["data"]]

# 使用 ggplot2 重绘精美气泡图
ggplot(rt, aes(id, features.plot)) +
  theme_bw() +
  coord_flip() +
  xlab("") +
  ylab(" ") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled), alpha = 0.6) +
  scale_colour_gradient(low = "navy", high = "#ff0000")

# 指定 KEGG 通路的 gmt 文件
gmtFile <- "kegg_pathway.gmt"
geneSet <- getGmt(gmtFile, geneIdType = SymbolIdentifier())

# --- Seurat v5 语法更新：使用 LayerData 提取 counts ---
exp <- as.data.frame(LayerData(pbmc1, assay = "RNA", layer = "counts"))

dimnames <- list(rownames(exp), colnames(exp))
mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

# --- 运行 ssGSEA 算法（适配新版 GSVA 1.50+ 语法） ---

# 1. 首先创建一个 ssGSEA 专用的参数对象
# 注意：新版中 abs.ranking 和 min.sz 等参数都在这个函数里设置
ssgsea_param <- ssgseaParam(
  exprData = mat, 
  geneSets = geneSet, 
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25, # ssGSEA 的默认权重参数
  normalize = TRUE
)

# 2. 调用 gsva 函数，并传入参数对象和并行计算设置
ssgseaScore <- gsva(
  ssgsea_param, 
  BPPARAM = MulticoreParam(workers = 12) # 开启 8 核心并行计算
)

# 定义归一化函数
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# 对 ssGSEA 结果进行归一化处理
ssgseaOut <- normalize(ssgseaScore)

# 合并评分信息与 ssGSEA 得分
rt <- merge(score, t(ssgseaOut), by = 0)
rt$score <- as.numeric(rt$score)

# 设置分组因子顺序
group <- levels(factor(rt$Type))
rt$Type <- factor(rt$Type, levels = group)

# 自动生成两两比较组合
comp <- combn(group, 2)
my_comparisons <- list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}

# 绘制最终的箱线图
ggplot(data = rt, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.5, size = 1.3) +
  geom_boxplot(alpha = 1, outlier.shape = NA, size = 1.3, width = 0.5, fatten = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(
    aes(group = Type),
    method = "wilcox.test",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
    label = "p.signif"
  ) +
  theme(legend.position = "right") +
  ggtitle("") +
  xlab("") +
  ylab("ICD score") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 40, hjust = 1)
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold", color = "black"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )
