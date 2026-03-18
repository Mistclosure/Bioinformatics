#C----
library(CRDscore)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggpubr) # 确保加载此包以支持 stat_compare_means
load("scRNA.Rdata")

# Seurat v5 语法：合并层级，确保后续计算一致性
scRNA <- JoinLayers(scRNA,assay = 'RNA')
# 1. 加载数据
load("Tcell.Rdata")
Tcell <- JoinLayers(Tcell,assay = 'RNA')
# --- 修改开始：替换 meat.txt 读取为标记法 ---
# 读取恶性细胞名单
exp <- as.data.frame(LayerData(scRNA, assay = "RNA", layer = "counts"))
rloop_data <- read.csv("Rloop.csv", header = T, check.names = F)
selected_genes <- rloop_data$GeneSymbol[1:92]
score <- cal_CRDscore(expr = exp, n.bins = 50, circadians = selected_genes, study.type = "scRNAseq")
score <- as.data.frame(score)
scRNA$score <- score$score
save(scRNA,file = 'scRNA.Rdata')
write.table(score, file = "CRDscore_all.txt", 
            sep = "\t", quote = FALSE, row.names = T)
idx <- match(colnames(Tcell), colnames(scRNA))
Tcell$score <- scRNA$score[idx]
# 在 Tcell 对象中直接标记 Type (Seurat v5 赋值语法)
Tcell$Type <- ifelse(Tcell$score > median(Tcell$score), "High", "Low")

# 提取元数据用于后续绘图
meta.data <- Tcell@meta.data
# --- 修改结束 ---

# 下面代码保持不动（计算比例与绘图）
a <- data.frame(table(meta.data$Type, meta.data$singleR))
a <- a[a$Var2!='Unclassified_T',]
a <- ddply(a, .(Var1), transform, percent = Freq/sum(Freq) * 100)
a$label = paste0(sprintf("%.1f", a$percent), "%")

# 定义一个包含 14 个颜色的向量，防止以后亚群变动再次报错
my_colors <- c("#8DD3C7", "#FFFFB3", "#FB8072", "#BEBADA", "#80B1D3", "#FDB462", 
               "#B3DE69", "#FCCDE5", "#BC80BD", "#D9D9D9", "#CCEBC5", "#FFED6F",
               "#7FC97F", "#BEAED4")

a %>%
  drop_na() %>%
  ggplot(aes(fill = Var2, y = percent, x = Var1)) +
  geom_bar(position = "fill", stat = "identity",width = 0.2) +
  # 使用我们定义的 14 色板
  scale_fill_manual(values = my_colors) + 
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent(%)") + xlab("") + labs(fill = "T cell") +
  theme_gray() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  )

# --- 第二部分：关联 Rloop 评分 ---
# 匹配细胞类型
scRNA$singleR = celltype[match(clusters, celltype$ClusterID), 'celltype']
# 1. 确保 singleR 是字符型，防止因 Factor 导致无法写入新标签
scRNA$singleR <- as.character(scRNA$singleR)

# 2. 将 Tcell 对象中细分过的标签按 Barcode 覆盖回 scRNA
# 这一步会把 scRNA 里的 "T cell" 变成具体的 "CD8 Exhaust" 等
scRNA@meta.data[colnames(Tcell), "singleR"] <- as.character(Tcell$singleR)
A = scRNA@meta.data
# --- 新增代码：按 orig.ident 分组计算 score 平均值 ---
# 1. 计算每个样本的平均分
sample_avg <- aggregate(score ~ orig.ident, data = A, FUN = mean)
colnames(sample_avg) <- c("SampleID", "Score")
write.table(sample_avg, file = "Rloopscore.txt", 
            sep = "\t", quote = FALSE, row.names = T)
# 2. 根据样本平均分的中位数划分为 High/Low (对接你图中的 Type)
sample_avg$Type <- ifelse(sample_avg$Score > median(sample_avg$Score), "High", "Low")
rownames(sample_avg) <- sample_avg$SampleID

rt = merge(sample_avg, data, by.x = 0, by.y = 0)
rt3 = rt[, c(4, 6)]
colnames(rt3) = c("Type", "Expression")
rt3$Expression = rt3$Expression * 100

group = levels(factor(rt3$Type))
rt3$Type = factor(rt3$Type, levels = group)
comp = combn(group, 2)
my_comparisons = list()
for (i in 1:ncol(comp)) { my_comparisons[[i]] <- comp[, i] }

# 绘制小提琴图
p = ggplot(rt3, aes(x = Type, y = Expression, fill = Type)) +
  geom_violin(trim = FALSE, color = "white") +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  geom_point(aes(x = Type, y = Expression), pch = 19, position = position_dodge(0.9), size = 0.5) +
  scale_fill_manual(values = c("#eaeae8", "#7bacb0")) +
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("") +
  ylab("T cell") + xlab("")

p

# 输出 PDF
pdf(file = "T cell.pdf", width = 2.5, height = 3.5)
plot(p)
dev.off()
