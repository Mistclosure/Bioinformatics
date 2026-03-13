#H----
library(GSVA)
library(limma)
library(GSEABase)
library(msigdbr)  # 新增：用于在线获取基因集
library(dplyr)
library(Seurat)

# 1. 加载 KEGG 基因集
m_df <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)
geneSet = split(m_df$gene_symbol, m_df$gs_name)

# 2. Seurat v5 语法提取数据
pbmc1 <- NormalizeData(pbmc1, assay = "RNA") 
exp = as.data.frame(LayerData(pbmc1, assay = "RNA", layer = "data"))

dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

# 3. 核心修改：适配新版 GSVA 语法，解决继承方法报错
# 创建参数对象 gsvaParam (针对原代码中的 method='gsva')
# 注意：参数名在新版中略有变化，如 abs.ranking 变为 absRanking, min.sz 变为 minSize

# 1. 正常的参数对象（这里不放并行参数）
param <- gsvaParam(exprData = mat, 
                   geneSets = geneSet, 
                   kcdf = "Gaussian", 
                   absRanking = TRUE, 
                   minSize = 10,
                   maxSize = 500) # 建议限制基因集最大值，也能提速

# 2. 引入并行计算包
library(BiocParallel)

# 3. 在 gsva() 函数中开启加速
# 设置 workers 为你想使用的 CPU 核心数
ssgseaScore = gsva(param, 
                   BPPARAM = MulticoreParam(workers = 32, progressbar = TRUE))

# --- 以下代码保持原样 ---
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)

# 4. 获取 Metadata
meta = pbmc1@meta.data
# 用这种方式代替 arrange，行名不会丢失
meta <- meta[order(meta$group, decreasing = TRUE), ]

GSVA_hall = ssgseaOut[,rownames(meta)]
GSVA_hall = GSVA_hall[-1,]
dimnames=list(rownames(GSVA_hall),colnames(GSVA_hall))
GSVA_hall=matrix(as.numeric(as.matrix(GSVA_hall)),nrow=nrow(GSVA_hall),
                 dimnames=dimnames)

# 5. Limma 差异分析
library(limma)
group <- factor(meta$group, levels = c('High', 'Low')) 
design <- model.matrix(~0+group)
colnames(design) = levels(group)
rownames(design) = colnames(GSVA_hall)

compare <- makeContrasts(High - Low, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)

# 修改：虽然提取 200 个，但后续我们会筛选最显著的来画图
dat_plot <- topTable(fit3, coef=1, number=200)

# 【核心修改：清洗通路名称】
# 去掉 REACTOME_ 前缀，并将下划线转为空格，首字母大写
dat_plot$id <- rownames(dat_plot)
dat_plot$id <- str_replace(dat_plot$id, "REACTOME_", "")
dat_plot$id <- str_replace_all(dat_plot$id, "_", " ")

# 6. 可视化准备
dat_plot$threshold = factor(ifelse(dat_plot$t >-2, ifelse(dat_plot$t >= 2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))

# --- 核心修改：清洗通路名称 (转为 Sentence case) ---
library(stringr)
dat_plot$id <- rownames(dat_plot) %>% 
  str_replace_all("KEGG_|REACTOME_|HALLMARK_", "") %>% # 去除前缀
  str_replace_all("_", " ") %>%                      # 下划线换空格
  str_to_sentence()                                 # 句首字母大写

# 确保 dat_plot 是标准 data.frame 并进行切片
dat_plot <- as.data.frame(dat_plot) %>% 
  dplyr::arrange(t) %>%
  dplyr::slice(c(1:15, (n()-14):n())) 

dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

library(ggplot2)
library(ggprism)

# 保存为 PDF
#pdf("GSVA_Reactome_T-value.pdf", width = 10, height = 8)

p <- ggplot(data = dat_plot, aes(x = id, y = t, fill = threshold)) +
  geom_col(width = 0.8) + 
  coord_flip() +
  # --- 核心修改：更换为图片中的配色 ---
  scale_fill_manual(values = c('Up' = '#f8766d', 'NoSignifi' = '#cccccc', 'Down' = '#66c2a5')) +
  geom_hline(yintercept = c(-2, 2), color = 'white', linewidth = 0.5, lty='dashed') +
  xlab('') +
  ylab('t value of GSVA score (High vs Low)') +
  guides(fill = "none") +
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(1,1,1,1), "cm") 
  )

# 动态获取当前展示的数量标签索引
low1 <- sum(dat_plot$t < 0)
high1 <- nrow(dat_plot)

# --- 核心修改：标签文字加粗 (fontface = "bold") ---
p + geom_text(data = dat_plot[1:low1,], aes(x = id, y = 0.5, label = id),
              hjust = 0, color = 'black', size = 3.5, fontface = "bold") + 
  geom_text(data = dat_plot[(low1 + 1):high1,], aes(x = id, y = -0.5, label = id),
            hjust = 1, color = 'black', size = 3.5, fontface = "bold")

dev.off()
