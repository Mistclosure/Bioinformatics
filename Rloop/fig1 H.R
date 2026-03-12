#H----
library(GSVA)
library(limma)
library(GSEABase)
library(msigdbr)  # 新增：用于在线获取基因集
library(dplyr)
library(Seurat)

load("Malignant.Rdata")

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
param <- gsvaParam(exprData = mat, 
                   geneSets = geneSet, 
                   kcdf = "Gaussian", 
                   absRanking = TRUE, 
                   minSize = 10)

# 执行打分
ssgseaScore = gsva(param)

# --- 以下代码保持原样 ---
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)

# 4. 获取 Metadata
meta = pbmc1@meta.data
meta = arrange(meta,desc(group))

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
dat_plot <- topTable(fit3, coef=1, number=200)
dat_plot$id <- rownames(dat_plot)

# 6. 可视化准备
library(stringr)
dat_plot$threshold = factor(ifelse(dat_plot$t >-2, ifelse(dat_plot$t >= 
                                                            2 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

library(ggplot2)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up' = '#36638a','NoSignifi' = '#cccccc','Down' = '#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',linewidth = 0.5, lty='dashed') +
  xlab('') +
  ylab('t value of GSVA score') +
  guides(fill=F)+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# 7. 添加文字标签
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
high1 <- nrow(dat_plot)
p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'black') +
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black')
