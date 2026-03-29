# ==============================================================================
# 1. T 细胞亚群提取与重聚类
# ==============================================================================

# 加载全局单细胞数据
load("scRNA.Rdata")

# 提取元数据并根据之前的 SingleR 鉴定结果筛选 T 细胞
meta = scRNA@meta.data
T_cell_names = rownames(meta[meta$singleR1 == "T_cells",])
T_cell = scRNA[, T_cell_names]

# 初步查看 T 细胞分布
DimPlot(T_cell, label.size = 5, reduction = 'umap')

# 【Seurat v5 语法】提取 T 细胞的原始计数矩阵并重新创建对象以进行亚群分析
exp = as.data.frame(LayerData(T_cell, assay = "RNA", layer = "counts"))
pbmc1 <- CreateSeuratObject(counts = exp, project = "T cell")

# 标准预处理流程：标准化、找高变基因、缩放
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc1)
pbmc1 <- ScaleData(pbmc1, features = all.genes)

# 线性与非线性降维
pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1))
pbmc1 <- FindNeighbors(pbmc1, dims = 1:20)
pbmc1 <- FindClusters(pbmc1, resolution = 1)
pbmc1 <- RunTSNE(pbmc1, dims = 1:20)
pbmc1 <- RunUMAP(pbmc1, dims = 1:20)

# 可视化亚群分布
DimPlot(pbmc1, reduction = "umap", label = T)
DimPlot(pbmc1, reduction = "tsne", label = T)

# 定义不同 T 细胞亚群的 Marker 基因（耗竭、毒性、极化、记忆、调节等）
select_genes <- c("NKG7","IFNG","GZMK","GZMB","PRF1","GNLY","GZMA","IL2", 
                  "LAG3","LAYN","CTLA4","PDCD1", "CCR7","LEF1","SELL","TCF7", "FOXP3","IL2RA")

# 气泡图检查 Marker 表达情况
DotPlot(pbmc1, features = select_genes) + coord_flip()

# 读取外部定义的 T 细胞细分注释文件并映射
celltype = read.table("T_cell.txt", header=T, sep="\t", check.names=F)
# 【Seurat v5 语法】使用 $ 直接访问或赋值 metadata
clusters <- pbmc1$seurat_clusters
pbmc1$singleR1 = celltype[match(clusters, celltype$ClusterID), 'celltype']

# 最终可视化 T 细胞亚群注释结果
DimPlot(pbmc1, group.by = "singleR1", label.size = 5, reduction = 'umap')
DimPlot(pbmc1, reduction = "umap", label = T)
save(pbmc1, file = 'Tcell.Rdata')

# ==============================================================================
# 2. 准备 CellPhoneDB 输入文件（构建 Rloop 高低分组）
# ==============================================================================

# 读取带有 Rloop 分组信息的元数据
meat = read.table("pbmc1@meta.data.txt", header=T, sep="\t", check.names=F, row.names = 1)

# 使用 Rcpp 定义一个高效的稀疏矩阵转密集矩阵函数（应对大数据集）
library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp, NumericVector cp, NumericVector z, int nrows, int ncols){
  int k = z.size() ;
  IntegerMatrix mat(nrows, ncols);
  for (int i = 0; i < k; i++){
    mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}')

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1, mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x, nrows = mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# 【Seurat v5 语法】提取全局 RNA 计数矩阵
exp = as_matrix(LayerData(scRNA, assay = "RNA", layer = "counts"))

# 分别提取 T 细胞和巨噬细胞（Mac）的元数据与对应的表达矩阵，准备 CellPhoneDB 所需的 counts 和 meta
# 此处逻辑涉及：将肿瘤细胞（High/Low）与 免疫细胞（T/Mac）合并
load("Tcell.Rdata")
Tmeat = pbmc1@meta.data # 此处 pbmc1 已是 T 细胞子集
load("Mac.Rdata")
Macmeat = pbmc1@meta.data # 假设 Mac.Rdata 载入的对象也叫 pbmc1

# 构建 High 分组的 CellPhoneDB 输入
high_meta = meat[meat$group == "High", 8, drop=F] # 提取分组列
Tmeat_sub = Tmeat[, 7, drop=F]; colnames(Tmeat_sub) = "group"
Macmeat_sub = Macmeat[, 7, drop=F]; colnames(Macmeat_sub) = "group"

meat_high_all = rbind(high_meta, Tmeat_sub, Macmeat_sub)
meat_high_all = cbind(rownames(meat_high_all), meat_high_all)
same_h = intersect(rownames(meat_high_all), colnames(exp))
meat_high_final = meat_high_all[same_h,]

# 输出 High 组文件
write.table(meat_high_final, file="meat_High.txt", quote=F, sep="\t", row.names = F)
exp_high = exp[, rownames(meat_high_final)]
write.table(cbind(Gene=rownames(exp_high), exp_high), file="count_High.txt", quote=F, sep="\t", row.names = F)

# 构建 Low 分组的 CellPhoneDB 输入（逻辑同上）
low_meta = meat[meat$group == "Low", 8, drop=F]
meat_low_all = rbind(low_meta, Tmeat_sub, Macmeat_sub)
meat_low_all = cbind(rownames(meat_low_all), meat_low_all)
same_l = intersect(rownames(meat_low_all), colnames(exp))
meat_low_final = meat_low_all[same_l,]

# 输出 Low 组文件
write.table(meat_low_final, file="meat_Low.txt", quote=F, sep="\t", row.names = F)
exp_low = exp[, rownames(meat_low_final)]
write.table(cbind(Gene=rownames(exp_low), exp_low), file="count_Low.txt", quote=F, sep="\t", row.names = F)

# ==============================================================================
# 3. CellPhoneDB 结果可视化（趋化因子、协同刺激、协同抑制）
# ==============================================================================

# 定义感兴趣的相互作用对类别
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair, value = T)
costimulatory <- grep("CD86|CD80|CD48|TNF|CD2|CD40|CD27|CD28|CD58|ICOS|CD[1-9]", mymeans$interacting_pair, value = T)
coinhibitory <- grep("SIRP|CD47|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", mymeans$interacting_pair, value = T)

# 可视化函数逻辑：读取 High/Low 结果并合并绘图
# 以趋化因子为例：
process_cpdb <- function(path, group_tag, genes) {
  mypvals = read.table(paste0(path, "/pvalues.txt"), header=T, sep="\t", check.names=F, row.names=1)
  mymeans = read.table(paste0(path, "/means.txt"), header=T, sep="\t", check.names=F, row.names=1)
  
  # 筛选与转换数据
  mymeans %>% dplyr::filter(interacting_pair %in% genes) %>%
    dplyr::select(interacting_pair, starts_with(group_tag), ends_with(group_tag)) %>%
    reshape2::melt() -> mdf
  
  mypvals %>% dplyr::filter(interacting_pair %in% genes) %>%
    dplyr::select(interacting_pair, starts_with(group_tag), ends_with(group_tag)) %>%
    reshape2::melt() -> pdf
  
  res = merge(pdf, mdf, by = c("interacting_pair", "variable"))
  colnames(res) = c("interacting_pair", "CC", "pvals", "means")
  res$Type = paste0("Rloop_", group_tag)
  return(res)
}

# 分别处理 High 和 Low 的趋化因子数据并合并
A = process_cpdb("High", "High", chemokines)
B = process_cpdb("Low", "Low", chemokines)
E = rbind(A, B) %>% filter(pvals < 0.05)

# 气泡图绘制
ggplot(E, aes(CC, interacting_pair)) + 
  facet_grid(. ~ Type) +
  geom_point(aes(size = -log10(pvals + 0.0001), color = means), alpha = 0.6) +
  scale_colour_gradientn(colours = c('blue','cyan','white','orange','red')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
