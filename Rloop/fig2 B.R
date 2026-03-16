# 加载所有必要的 R 包
library(Seurat)       # 处理单细胞数据核心包
library(ggplot2)      # 绘图核心包
library(dplyr)        # 数据清洗与处理
library(magrittr)     # 提供管道操作符 %>%
library(RColorBrewer) # 提供 scale_fill_brewer 所需的调色板
library(CRDscore)
# B----
load("scRNA.Rdata")

# --- 修改部分：从 Rloop.csv 选取前 92 个基因 ---
rloop_data <- read.csv("Rloop.csv", header = T, check.names = F)
selected_genes <- rloop_data$GeneSymbol[1:92]
scRNA <- JoinLayers(scRNA, assay = "RNA")
# --- Seurat v5 语法更新：使用 LayerData 提取 counts ---
exp <- as.data.frame(LayerData(scRNA, assay = "RNA", layer = "counts"))

# 计算评分，使用选取的 92 个基因
score <- cal_CRDscore(expr = exp, n.bins = 50, circadians = selected_genes, study.type = "scRNAseq")
score <- as.data.frame(score)

A <- scRNA@meta.data
B <- as.data.frame(table(A$orig.ident))
score$id <- A$orig.ident
id <- c()
exp_val <- c() # 原变量名为 exp，与矩阵重名，建议保持逻辑但注意区分

for (i in B[, 1]) {
  rt <- score[score$id == i, ]
  sum_val <- log2(sum(rt$score))
  id <- c(id, i)
  exp_val <- c(exp_val, sum_val)
}

data <- data.frame(id, exp = exp_val)
data$Type <- ifelse(data$exp > median(data$exp), "High", "Low")
l <- data[data$Type == "Low", ]
h <- data[data$Type == "High", ]

Low <- B[which(B$orig.ident %in% l[, 1]), ]
Low$Type <- "Low"
High <- B[which(B$orig.ident %in% h[, 1]), ]
High$Type <- "High"

A1 <- rbind(High, Low)
A2 <- A1

# 细胞类型归类
A2$singleR1[which(A2$singleR1 == 'B_cell')] <- 'immune cell'
A2$singleR1[which(A2$singleR1 == 'DC')] <- 'immune cell'
A2$singleR1[which(A2$singleR1 == 'Endothelial_cells')] <- 'Other'
A2$singleR1[which(A2$singleR1 == 'Epithelial_cells')] <- 'Other'
A2$singleR1[which(A2$singleR1 == 'Fibroblasts')] <- 'Other'
A2$singleR1[which(A2$singleR1 == 'Macrophage')] <- 'immune cell'
A2$singleR1[which(A2$singleR1 == 'Mast')] <- 'immune cell'
A2$singleR1[which(A2$singleR1 == 'NK_cell')] <- 'immune cell'
A2$singleR1[which(A2$singleR1 == 'Oligodendrocytes')] <- 'Other'
A2$singleR1[which(A2$singleR1 == 'T_cells')] <- 'immune cell'

Low <- A2[A2$Type == "Low", ]
High <- A2[A2$Type == "High", ]

# 绘制 Low 组饼图
# 
diamonds_df <- Low %>%
  group_by(Type, singleR1) %>%
  tally() %>%
  mutate(prop = round((n / sum(n)) * 100, digits = 2)) %>%
  ungroup()

ggplot(diamonds_df, aes(x = Type, y = prop, fill = singleR1)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_brewer(palette = "Set3") +
  coord_polar("y", start = 0) +
  facet_wrap(~Type) +
  theme_void() +
  theme(strip.text = element_text(size = 16), legend.position = "top")

# 绘制 High 组饼图
diamonds_df <- High %>%
  group_by(Type, singleR1) %>%
  tally() %>%
  mutate(prop = round((n / sum(n)) * 100, digits = 2)) %>%
  ungroup()

ggplot(diamonds_df, aes(x = Type, y = prop, fill = singleR1)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_brewer(palette = "Set3") +
  coord_polar("y", start = 0) +
  facet_wrap(~Type) +
  theme_void() +
  theme(strip.text = element_text(size = 16), legend.position = "top")
