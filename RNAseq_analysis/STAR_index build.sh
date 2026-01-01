##这里以人类基因组为例，请按实际修改
cd /home/qiuzerui/annotationHv49

cat << 'EOF' > prepare_indices.sh
#!/bin/bash

# ==========================================
# 准备工作：解压与构建索引 (PE150优化版)
# ==========================================

# 设置路径
WORKDIR="/home/qiuzerui/annotationHv49"
STAR_OUT_DIR="/home/qiuzerui/star_index_h38"

# 文件名定义
FILE_FA_GZ="h38_PRI.fa.gz"
FILE_GENE_GTF_GZ="gencode.v49.annotation_PRI.gtf.gz"
FILE_TE_GZ="h38_TE.gz"

# 解压后的目标文件名
FASTA="h38_PRI.fa"
GENE_GTF="gencode.v49.annotation_PRI.gtf"
TE_GTF="h38_TE.gtf"

echo "=== 1. 解压文件 ==="

# 解压 FASTA
if [ ! -f "${FASTA}" ]; then
    if [ -f "${FILE_FA_GZ}" ]; then
        echo "正在解压基因组 FASTA..."
        gunzip -k "${FILE_FA_GZ}"
    else
        echo "错误：找不到 ${FILE_FA_GZ}"
        exit 1
    fi
else
    echo "${FASTA} 已存在，跳过解压。"
fi

# 解压 Gene GTF
if [ ! -f "${GENE_GTF}" ]; then
    if [ -f "${FILE_GENE_GTF_GZ}" ]; then
        echo "正在解压 Gene GTF..."
        gunzip -k "${FILE_GENE_GTF_GZ}"
    else
        echo "错误：找不到 ${FILE_GENE_GTF_GZ}"
        exit 1
    fi
else
    echo "${GENE_GTF} 已存在，跳过解压。"
fi

# 解压 TE GTF (重命名为 .gtf)
if [ ! -f "${TE_GTF}" ]; then
    if [ -f "${FILE_TE_GZ}" ]; then
        echo "正在解压 TE GTF..."
        gunzip -c "${FILE_TE_GZ}" > "${TE_GTF}"
    else
        echo "错误：找不到 ${FILE_TE_GZ}"
        exit 1
    fi
else
    echo "${TE_GTF} 已存在，跳过解压。"
fi


echo "=== 2. 构建 rRNA 和 mtDNA 索引 (Bowtie2) ==="
# 输出前缀: rRNA_mtDNA_index

# 2.1 提取 chrM
echo "提取 chrM (Mitochondria)..."
if ! command -v samtools &> /dev/null; then
    echo "错误: 未找到 samtools，请先安装 (conda install samtools)"
    exit 1
fi
samtools faidx "${FASTA}" chrM > chrM.fa

# 2.2 尝试提取 rRNA
echo "尝试提取 rRNA 序列..."
if command -v bedtools &> /dev/null; then
    grep 'gene_type "rRNA"' "${GENE_GTF}" | awk '$3=="exon"' > rRNA.temp.gtf
    if [ -s rRNA.temp.gtf ]; then
        bedtools getfasta -fi "${FASTA}" -bed rRNA.temp.gtf -fo rRNA.fa
        cat chrM.fa rRNA.fa > contamination.fa
        echo "成功提取 rRNA + chrM。"
    else
        echo "警告: GTF中未找到 'gene_type \"rRNA\"'，将仅使用 chrM。"
        cp chrM.fa contamination.fa
    fi
    rm -f rRNA.temp.gtf
else
    echo "警告: 未找到 bedtools，跳过 rRNA 提取，仅使用 chrM。"
    cp chrM.fa contamination.fa
fi

# 2.3 建立 Bowtie2 索引
echo "正在构建 Bowtie2 去污染索引..."
bowtie2-build --threads 12 contamination.fa rRNA_mtDNA_index


echo "=== 3. 构建 STAR 索引 (仅使用 Gene GTF) ==="
# 输出目录: /home/qiuzerui/star_index_h38

mkdir -p "${STAR_OUT_DIR}"

echo "使用基因组: ${FASTA}"
echo "使用注释: ${GENE_GTF}"
echo "优化参数: sjdbOverhang = 149 (针对 PE150)"

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir ${STAR_OUT_DIR} \
     --genomeFastaFiles ${FASTA} \
     --sjdbGTFfile ${GENE_GTF} \
     --sjdbOverhang 149 \
     --genomeSAsparseD 3

echo "=========================================="
echo "准备工作完成！"
echo "1. STAR 索引已生成于: ${STAR_OUT_DIR}"
echo "2. 去污染索引已生成于: ${WORKDIR}/rRNA_mtDNA_index"
echo "3. TEcount 所需 Gene GTF: ${WORKDIR}/${GENE_GTF}"
echo "4. TEcount 所需 TE GTF:   ${WORKDIR}/${TE_GTF}"
echo "=========================================="
EOF

chmod +x prepare_indices.sh
./prepare_indices.sh
