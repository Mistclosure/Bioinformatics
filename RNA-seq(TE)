#!/bin/bash

# ==========================================
# TE Analysis Pipeline (Mouse m39/GRCm39 完整版)
# ==========================================
# 功能：SRA转码 -> 质控 -> 去rRNA -> 比对 -> 定量

# 检查环境
for cmd in fastq-dump fastp bowtie2 STAR TEcount; do
    if ! command -v $cmd &> /dev/null; then
        echo "错误: 命令 $cmd 未找到！请先运行: conda activate te_env (或对应环境)"
        # 注意：如果没有 fasterq-dump，脚本会自动尝试 fastq-dump，但建议安装 sra-tools
    fi
done

# ======================
# 1. 路径配置
# ======================
# 工作主目录
WORKDIR="/home/qiuzerui/Phf8"
cd ${WORKDIR} || { echo "错误: 无法进入目录 ${WORKDIR}"; exit 1; }

# 输入目录
SRA_DIR="${WORKDIR}/sra"

# 输出目录结构
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 外部资源路径 (根据您的附件和要求修改)
STAR_INDEX="/home/qiuzerui/star_index_m39"
GTF_GENE="/home/qiuzerui/annotationMv38/gencode.vM38.annotation_PRI.gtf"
GTF_TE="/home/qiuzerui/annotationMv38/m39_TE.gtf"
RRNA_INDEX="/home/qiuzerui/annotationMv38/rRNA_mtDNA_index"

# 创建所有必要的输出目录
mkdir -p ${RAW_DIR} ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}

# ==========================================
# Step 0: SRA 转 FASTQ (新增步骤)
# ==========================================
echo "=== 开始 Step 0: SRA 格式转换 ==="

# 检查是否有 sra 文件
count_sra=$(ls ${SRA_DIR}/*.sra 2>/dev/null | wc -l)
if [ "$count_sra" -eq "0" ]; then
    echo "警告: 在 ${SRA_DIR} 未找到 .sra 文件。"
else
    for sra_file in ${SRA_DIR}/*.sra; do
        [ -e "$sra_file" ] || continue
        
        filename=$(basename ${sra_file})
        sample_name=${filename%.sra}
        
        echo ">>> 正在转换: ${sample_name} <<<"

        # 检查目标文件是否已存在 (避免重复转换)
        if [ -f "${RAW_DIR}/${sample_name}_1.fastq.gz" ]; then
            echo "   跳过 (FASTQ 已存在)"
            continue
        fi

        # 优先使用 fasterq-dump (多线程，速度快)，如果不行则回退到 fastq-dump
        if command -v fasterq-dump &> /dev/null; then
            # --split-3 确保双端测序被正确拆分
            # -e 8 使用8个线程
            fasterq-dump --split-3 --threads 8 --outdir ${RAW_DIR} --progress ${sra_file}
        else
            echo "   未找到 fasterq-dump，使用较慢的 fastq-dump..."
            fastq-dump --split-3 --gzip --outdir ${RAW_DIR} ${sra_file}
        fi

        # 如果使用的是 fasterq-dump，产生的是未压缩的 .fastq，需要手动 gzip
        if [ -f "${RAW_DIR}/${sample_name}_1.fastq" ]; then
            echo "   正在压缩 FASTQ 文件..."
            gzip -f ${RAW_DIR}/${sample_name}_1.fastq
            gzip -f ${RAW_DIR}/${sample_name}_2.fastq
        fi
    done
fi


# ==========================================
# Step 1-3: 预处理和比对循环
# ==========================================
echo "=== 开始 Step 1-3: 预处理和比对 ==="

# 注意：SRA工具生成的标准后缀通常是 _1.fastq.gz 和 _2.fastq.gz
# 之前的代码是 _f1 / _r2，这里统一修改为适配 SRA 工具的输出
for r1_file in ${RAW_DIR}/*_1.fastq.gz
do
    [ -e "$r1_file" ] || continue

    filename=$(basename ${r1_file})
    # 截取样本名 (去除 _1.fastq.gz)
    sample_name=${filename%_1.fastq.gz}
    r2_file="${RAW_DIR}/${sample_name}_2.fastq.gz"

    echo "------------------------------------------"
    echo ">>> 正在处理: ${sample_name} <<<"
    echo "------------------------------------------"

    # Step 1: Fastp
    echo "[1/3] 运行 fastp..."
    if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
        if [ ! -f "$r2_file" ]; then
             echo "错误: 找不到对应的 R2 文件: $r2_file"
             continue
        fi

        fastp -i ${r1_file} -I ${r2_file} \
              -o ${TRIM_DIR}/${sample_name}_1.clean.fq.gz \
              -O ${TRIM_DIR}/${sample_name}_2.clean.fq.gz \
              -h ${TRIM_DIR}/${sample_name}_fastp.html \
              -j ${TRIM_DIR}/${sample_name}_fastp.json \
              --thread 8 --detect_adapter_for_pe --length_required 25
    else
        echo "   跳过 (文件已存在)"
    fi

    
    # Step 2: Bowtie2 (去rRNA/mtDNA)
    echo "[2/3] 去除 rRNA 和 线粒体..."
    if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
        # 检查 rRNA 索引完整性 (检查其中一个后缀即可)
        if [ ! -f "${RRNA_INDEX}.1.bt2" ] && [ ! -f "${RRNA_INDEX}.1.bt2l" ]; then
             echo "错误: rRNA索引未找到: ${RRNA_INDEX}"
             exit 1
        fi

        bowtie2 -p 12 --very-fast-local --no-unal -x ${RRNA_INDEX} \
                -1 ${TRIM_DIR}/${sample_name}_1.clean.fq.gz \
                -2 ${TRIM_DIR}/${sample_name}_2.clean.fq.gz \
                --un-conc-gz ${CLEAN_DIR}/${sample_name}_clean \
                > /dev/null
        
        # Bowtie2输出的 _clean.1.gz 重命名为 _1.final.fq.gz
        # 注意：bowtie2 输出压缩文件时后缀可能是 .gz
        mv ${CLEAN_DIR}/${sample_name}_clean.1 ${CLEAN_DIR}/${sample_name}_1.final.fq.gz 2>/dev/null || mv ${CLEAN_DIR}/${sample_name}_clean.1.gz ${CLEAN_DIR}/${sample_name}_1.final.fq.gz
        mv ${CLEAN_DIR}/${sample_name}_clean.2 ${CLEAN_DIR}/${sample_name}_2.final.fq.gz 2>/dev/null || mv ${CLEAN_DIR}/${sample_name}_clean.2.gz ${CLEAN_DIR}/${sample_name}_2.final.fq.gz
    else
        echo "   跳过 (文件已存在)"
    fi

    # Step 3: STAR
    echo "[3/3] 运行 STAR..."
    if [ ! -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
        STAR --runThreadN 12 --genomeDir ${STAR_INDEX} \
             --readFilesIn ${CLEAN_DIR}/${sample_name}_1.final.fq.gz ${CLEAN_DIR}/${sample_name}_2.final.fq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix ${ALIGN_DIR}/${sample_name}. \
             --outSAMtype BAM SortedByCoordinate \
             --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
             --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All

        samtools index ${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam
    else
        echo "   跳过 (BAM 已存在)"
    fi
done

# ==========================================
# Step 4: TEcount
# ==========================================
echo "=== 开始 Step 4: TEcount 定量 ==="

for bam_file in ${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam
do
    [ -e "$bam_file" ] || continue
    sample_name=$(basename "${bam_file}" .Aligned.sortedByCoord.out.bam)
    
    echo "正在定量样本: ${sample_name}"
    
    if [ ! -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then
        TEcount --sortByPos --format BAM --mode multi \
                --GTF "${GTF_GENE}" \
                --TE "${GTF_TE}" \
                --project "${COUNTS_DIR}/${sample_name}" \
                --stranded reverse \
                -b "${bam_file}"
    else
        echo "-> 跳过 (结果已存在)"
    fi
done
echo "所有流程完成。"
