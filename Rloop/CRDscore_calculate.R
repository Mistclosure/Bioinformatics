library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(tidyverse)
library(plyr)
library(Seurat)
library(CRDscore)
setwd('/mnt/disk1/qiuzerui/downloads/Rloop/')
load("Malignant.Rdata")
pbmc1 <- JoinLayers(pbmc1,assay = 'RNA')
pbmc1 <- NormalizeData(pbmc1,assay = 'RNA')
# 1. Seurat v5 提取表达矩阵
# 如果你的对象是合并过的，v5 建议先 JoinLayers 确保数据完整

# 使用 LayerData 提取 counts (等同于 v4 的 @assays$RNA@counts)
exp = as.data.frame(LayerData(pbmc1, assay = "RNA", layer = "data"))

# 2. 提取 Rloop.csv 中的基因集
rloop_data = read.csv("genes2.csv", header = T, check.names = F)
# 提取 GeneSymbol 列的前 92 个基因
target_genes = as.character(rloop_data$Gene[1:100])
# 建议加上交集判断，防止打分函数因为找不到基因而报错
target_genes = intersect(target_genes, rownames(pbmc1))

# 3. 计算评分
score <- cal_CRDscore(expr = exp, n.bins = 50, circadians = 
                        target_genes, study.type = "scRNAseq")

score = as.data.frame(score)
score = cbind(id=rownames(score), score)
# 核心修改：score 的拆分处理（格式与 meta 一致）
score$id = rownames(score)
score$Original_Sample_ID = sub("^[^_]*_", "", rownames(score))

# 提取 metadata 并添加细胞 ID 列（id）和拆分出的患者 ID 列（Original_Sample_ID）
meta = pbmc1@meta.data
meta$id = rownames(meta) 
meta$Original_Sample_ID = sub("^[^_]*_", "", rownames(meta))

# 读取 Cli.csv 文件
cli = read.csv("Cli.csv", header=T, check.names=F)

# 1. 按照患者 ID 合并细胞元数据与临床信息
rt = merge(meta, cli, by.x= "Type",by.y="ID")

# 2. 按照细胞 ID 合并评分结果与上一步的合并表
# 合并依旧按照指定列名的方式
rt1 = merge(score, rt, by.x = "id", by.y = "id")

# --- 绘图部分 1 (Column 12) ---
# 注意：由于 merge 后列顺序改变，建议检查 rt1 的列名以确保索引正确
data = rt1[,c(2,17)]
colnames(data) = c("score","Type")
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=group)

# 增加安全检查防止 combn 报错
my_comparisons=list()
if(length(group) >= 2){
  comp=combn(group,2)
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
}

ggplot(data,aes(x = Type , y=score, color=Type))+
  stat_boxplot(geom="errorbar",width=0.6)+
  geom_boxplot(alpha=0.7,outlier.shape = NA,size=0.7,width=0.7,fatten=0.7)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
  theme(legend.position = "none")+
  ggtitle("") +
  ylim(-0.4, 0.4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("RloopScore")+
  theme(axis.title.x=element_blank(),
        axis.text.y = element_text(size = 15,face = "bold",color = "black"),
        axis.title.y = element_text(size = 15,face = "bold",color = "black"),
        axis.text.x = element_text(size = 15,face = "bold",color = "black"),
        legend.position = "none",
        plot.title = element_text(face = "bold",size=15,hjust = 0.5))
