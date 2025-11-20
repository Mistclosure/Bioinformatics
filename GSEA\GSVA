# ==============================================================================
# 最终版 (专业绘图优化)：TE/基因 通路分析流程
# 1. ID 转换：使用 biomaRt (适配 Gencode v49)
# 2. GSEA图：P值显著性配色 (红=显著)
# 3. GSVA图：
#    - 开启行标准化 (Z-score)
#    - 【优化】样本标签改为具体的 "sgNC", "sgKAT8"
#    - 【修复】分组颜色显式指定 (Treat=红, Control=蓝)
# 4. 稳定性：强制单核运行
# ==============================================================================

# 1. 环境准备
# ------------------------------------------------------------------------------
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(GSVA)
library(GSEABase)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(pheatmap)
library(msigdbr)
library(BiocParallel)

if (!require("biomaRt", quietly = TRUE)) {
  stop("请先安装 biomaRt: BiocManager::install('biomaRt')")
}

register(SerialParam())

# ==============================================================================
# 2. 数据读取与清洗
# ==============================================================================
counts_df <- read.csv("KAT8_counts.csv", row.names = 1)

is_gene <- str_detect(rownames(counts_df), "^ENSG")
gene_counts <- counts_df[is_gene, ]

cleaned_ids <- str_split(rownames(gene_counts), "\\.", simplify = TRUE)[, 1]
rownames(gene_counts) <- cleaned_ids

# ==============================================================================
# 2.1 ID 转换 (BiomaRt)
# ==============================================================================
print("--- 开始 ID 转换 ---")
gene_map <- NULL

tryCatch({
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
  print("正在查询 Ensembl 数据库...")
  raw_map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"),
                   filters = "ensembl_gene_id",
                   values = rownames(gene_counts),
                   mart = mart)
  
  raw_map <- raw_map[!is.na(raw_map$entrezgene_id), ]
  coding_genes <- raw_map[raw_map$gene_biotype == "protein_coding", ]
  
  gene_map <- data.frame(
    ENSEMBL = coding_genes$ensembl_gene_id,
    ENTREZID = as.character(coding_genes$entrezgene_id),
    SYMBOL = coding_genes$hgnc_symbol,
    stringsAsFactors = FALSE
  )
  print(paste("BiomaRt 成功! 保留蛋白编码基因:", length(unique(coding_genes$ensembl_gene_id))))
  
}, error = function(e) {
  stop(paste("BiomaRt 连接失败:", e$message))
})

gene_counts <- gene_counts[rownames(gene_counts) %in% gene_map$ENSEMBL, ]
gene_counts$ENSEMBL <- rownames(gene_counts)
merged_data <- merge(gene_counts, gene_map, by = "ENSEMBL")
merged_data <- merged_data[!duplicated(merged_data$ENTREZID), ]
rownames(merged_data) <- merged_data$ENTREZID

expr_matrix <- as.matrix(merged_data[, 2:5]) 
colnames(expr_matrix) <- colnames(counts_df)[1:4]

# ==============================================================================
# 3. 差异分析
# ==============================================================================
# 定义样本信息
# 91,92 -> sgNC (Control)
# 93,94 -> sgKAT8 (Treat)
sample_info <- data.frame(
  row.names = colnames(expr_matrix),
  condition = factor(c("Control", "Control", "Treat", "Treat"), 
                     levels = c("Control", "Treat"))
)

dds <- DESeqDataSetFromMatrix(countData = expr_matrix, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

# 明确对比：Treat (sgKAT8) vs Control (sgNC)
# Log2FC > 0 表示 sgKAT8 上调
res <- results(dds, contrast = c("condition", "Treat", "Control"))

gene_list <- res$stat
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

# ==============================================================================
# 4. 关注的通路定义
# ==============================================================================
target_kegg_ids <- c("hsa04622", "hsa04620", "hsa04623", "hsa04621", "hsa04612")
target_go_ids   <- c("GO:0051607", "GO:0060337")

kegg_id_map <- c(
  "hsa04622" = "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "hsa04620" = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "hsa04623" = "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
  "hsa04621" = "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "hsa04612" = "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

go_id_map <- c(
  "GO:0051607" = "GOBP_DEFENSE_RESPONSE_TO_VIRUS",
  "GO:0060337" = "GOBP_TYPE_I_INTERFERON_SIGNALING_PATHWAY"
)

target_ids_all <- c(as.character(kegg_id_map), target_go_ids)

pathway_names <- c(
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY" = "RIG-I-like receptor signaling",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"  = "Toll-like receptor signaling",
  "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"         = "Cytosolic DNA-sensing",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"   = "NOD-like receptor signaling",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"   = "Antigen processing & presentation",
  "GO:0051607" = "Defense response to virus",
  "GO:0060337" = "Type I interferon signaling",
  "GOBP_DEFENSE_RESPONSE_TO_VIRUS"             = "Defense response to virus",
  "GOBP_TYPE_I_INTERFERON_SIGNALING_PATHWAY"   = "Type I interferon signaling"
)

# ==============================================================================
# 5. GSEA 分析
# ==============================================================================
print("--- 开始 GSEA 分析 ---")
m_df_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
m_df_kegg <- m_df_c2 %>% dplyr::filter(grepl("^KEGG_", gs_name))
kegg_t2g <- m_df_kegg %>% dplyr::select(gs_name, entrez_gene)
kk_gsea <- GSEA(gene_list, TERM2GENE = kegg_t2g, pvalueCutoff = 1, verbose = FALSE, seed = 123)

go_gsea <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", 
                 pvalueCutoff = 1, minGSSize = 5, verbose = FALSE, seed = 123)

res_df_kegg <- as.data.frame(kk_gsea)
res_df_go   <- as.data.frame(go_gsea)
all_gsea_res <- rbind(res_df_kegg, res_df_go)

print("正在导出 GSEA 原始结果 (GSEA_All_Results_Unfiltered.csv)...")
write.csv(all_gsea_res, "GSEA_All_Results_Unfiltered.csv", row.names = FALSE)

plot_data <- all_gsea_res[all_gsea_res$ID %in% target_ids_all, ]

if(nrow(plot_data) > 0) {
  plot_data$Description_Label <- pathway_names[plot_data$ID]
  plot_data$Description_Label <- ifelse(is.na(plot_data$Description_Label), plot_data$ID, plot_data$Description_Label)
  plot_data <- plot_data[order(plot_data$NES), ]
  plot_data$Description_Label <- factor(plot_data$Description_Label, levels = plot_data$Description_Label)
  
  p_dot <- ggplot(plot_data, aes(x = NES, y = Description_Label)) +
    geom_segment(aes(x = 0, xend = NES, y = Description_Label, yend = Description_Label), color = "grey70", size = 0.8) +
    geom_point(aes(color = p.adjust, size = setSize)) + 
    scale_color_gradient(low = "red", high = "blue") + # 显著红，不显著蓝
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(title = "Target Pathways Enrichment (sgKAT8 vs sgNC)", x = "NES", y = NULL)
  
  ggsave("GSEA_Summary_Dotplot_Final.pdf", p_dot, width = 8, height = 0.5 * nrow(plot_data) + 3)
  print("GSEA 气泡图已保存。")
}

# ==============================================================================
# 6. GSVA 分析 (专业绘图版)
# ==============================================================================
print("--- 开始 GSVA 分析 ---")
vst_data <- vst(dds, blind = FALSE)
expr_normalized <- assay(vst_data)

custom_gene_sets <- list()

target_kegg_names <- as.character(kegg_id_map)
for (kname in target_kegg_names) {
  genes <- m_df_kegg %>% filter(gs_name == kname) %>% pull(entrez_gene) %>% unique()
  if(length(genes) > 0) custom_gene_sets[[kname]] <- genes
}

m_df_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
for (gid in target_go_ids) {
  genes <- NULL
  key_used <- gid
  
  genes <- m_df_c5 %>% filter(gs_exact_source == gid) %>% pull(entrez_gene) %>% unique()
  if(length(genes) == 0) {
    target_name <- go_id_map[gid]
    if(!is.na(target_name)) {
      genes <- m_df_c5 %>% filter(gs_name == target_name) %>% pull(entrez_gene) %>% unique()
      if(length(genes) > 0) key_used <- target_name
    }
  }
  if(length(genes) == 0) {
    try({
      genes_org <- select(org.Hs.eg.db, keys = gid, keytype = "GOALL", columns = "ENTREZID")$ENTREZID
      genes <- unique(genes_org[!is.na(genes_org)])
      if(length(genes) > 0) key_used <- gid
    }, silent = TRUE)
  }
  if(length(genes) > 0) custom_gene_sets[[key_used]] <- genes
}

if(length(custom_gene_sets) > 0) {
  
  if (packageVersion("GSVA") >= "1.52.0") {
    gsva_par <- gsvaParam(exprData = expr_normalized, geneSets = custom_gene_sets, kcdf = "Gaussian", minSize = 1)
    gsva_res <- gsva(gsva_par, BPPARAM = SerialParam())
  } else {
    gsva_res <- gsva(expr_normalized, custom_gene_sets, kcdf = "Gaussian", min.sz = 1, parallel.sz = 1)
  }
  
  if(is.list(gsva_res) && !is.matrix(gsva_res) && "es" %in% names(gsva_res)) gsva_res <- gsva_res$es
  
  rn <- rownames(gsva_res)
  new_rn <- pathway_names[rn]
  rownames(gsva_res) <- ifelse(is.na(new_rn), rn, new_rn)
  
  print("正在导出 GSVA 评分矩阵 (GSVA_Scores_Matrix.csv)...")
  write.csv(as.data.frame(gsva_res), "GSVA_Scores_Matrix.csv", row.names = TRUE)
  
  save_height <- nrow(gsva_res) * 0.5 + 2
  save_width  <- ncol(gsva_res) * 0.5 + 5
  
  while(!is.null(dev.list())) dev.off()
  
  tryCatch({
    # === 专业绘图设置 ===
    
    # 1. 自定义分组颜色：Treat=红, Control=蓝
    ann_colors <- list(
      condition = c(Control = "blue", Treat = "red")
    )
    
    # 2. 【新增】自定义样本标签：使用真实的实验名称
    new_labels <- c("sgNC_1", "sgNC_2", "sgKAT8_1", "sgKAT8_2")
    
    # 3. 颜色标尺
    colorPalette <- colorRampPalette(c("blue", "white", "red"))(100)
    
    pheatmap(gsva_res, 
             cluster_cols = FALSE, cluster_rows = TRUE,
             main = "GSVA Scores (sgKAT8 vs sgNC)",
             annotation_col = sample_info[, "condition", drop=FALSE],
             
             annotation_colors = ann_colors,
             labels_col = new_labels, # 应用新标签
             
             cellheight = 20, cellwidth = 30, fontsize = 12,
             border_color = "white",
             scale = "row", 
             breaks = seq(-2, 2, length.out = 101),
             color = colorPalette,
             
             filename = "GSVA_Heatmap_Final.pdf",
             width = save_width, height = save_height)
    print("GSVA 热图已保存 (已应用 sgNC/sgKAT8 标签)。")
  }, error = function(e) { print(paste("绘图失败:", e$message)) })
  
} else {
  print("未构建出有效的基因集，GSVA 跳过。")
}
