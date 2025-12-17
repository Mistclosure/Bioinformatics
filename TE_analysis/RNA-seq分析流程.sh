#!/bin/bash

# ==========================================
# TE Analysis Pipeline (智能后缀识别版)
# ==========================================
# 功能：自动检测SRA -> (可选转换) -> 兼容fastq.gz/fq.gz -> 质控 -> 去rRNA -> 比对 -> 定量

# 检查环境
for cmd in fastq-dump fastp bowtie2 STAR TEcount; do
    if ! command -v $cmd &> /dev/null; then
        echo "错误: 命令 $cmd 未找到！请先运行: conda activate te_env (或对应环境)"
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

# 外部资源路径
STAR_INDEX="/home/qiuzerui/star_index_m39"
GTF_GENE="/home/qiuzerui/annotationMv38/gencode.vM38.annotation_PRI.gtf"
GTF_TE="/home/qiuzerui/annotationMv38/m39_TE.gtf"
RRNA_INDEX="/home/qiuzerui/annotationMv38/rRNA_mtDNA_index"

# 创建所有必要的输出目录
mkdir -p ${RAW_DIR} ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}

# ==========================================
# Step 0: SRA 转 FASTQ (智能检测)
# ==========================================
echo "=== Step 0: 检查 SRA 数据源 ==="

DO_CONVERSION=false

if [ -d "${SRA_DIR}" ]; then
    count_sra=$(ls ${SRA_DIR}/*.sra 2>/dev/null | wc -l)
    if [ "$count_sra" -gt "0" ]; then
        echo "✅ 检测到 ${SRA_DIR} 中有 $count_sra 个 SRA 文件。"
        DO_CONVERSION=true
    else
        echo "⚠️  检测到 ${SRA_DIR} 文件夹，但里面没有 .sra 文件。"
    fi
else
    echo "⚠️  未检测到 SRA 文件夹 (${SRA_DIR})。"
fi

if [ "$DO_CONVERSION" = true ]; then
    echo ">>> 开始 SRA -> FASTQ 转换..."
    for sra_file in ${SRA_DIR}/*.sra; do
        [ -e "$sra_file" ] || continue
        
        filename=$(basename ${sra_file})
        sample_name=${filename%.sra}
        
        echo "   正在处理: ${sample_name}"

        # 检查目标文件是否已存在 (检查 fastq.gz)
        if [ -f "${RAW_DIR}/${sample_name}_1.fastq.gz" ] || [ -f "${RAW_DIR}/${sample_name}_1.fq.gz" ]; then
            echo "   -> 跳过 (FASTQ/FQ 已存在)"
            continue
        fi

        # 转换逻辑
        if command -v fasterq-dump &> /dev/null; then
            fasterq-dump --split-3 --threads 8 --outdir ${RAW_DIR} --progress ${sra_file}
        else
            echo "   -> 使用 fastq-dump (较慢)..."
            fastq-dump --split-3 --gzip --outdir ${RAW_DIR} ${sra_file}
        fi

        # 压缩并标准化后缀 (默认统一压缩为 .fastq.gz，方便后续处理)
        if [ -f "${RAW_DIR}/${sample_name}_1.fastq" ]; then
            echo "   -> 正在压缩 FASTQ..."
            gzip -f ${RAW_DIR}/${sample_name}_1.fastq
            gzip -f ${RAW_DIR}/${sample_name}_2.fastq
        fi
    done
else
    echo "⏭️  跳过 Step 0。假设原始 FASTQ/FQ 文件已存在于: ${RAW_DIR}"
fi


# ==========================================
# Step 1-3: 预处理和比对循环 (支持 fastq.gz 和 fq.gz)
# ==========================================
echo "=== 开始 Step 1-3: 预处理和比对 ==="

# 启用 nullglob 选项，防止没有匹配文件时 for 循环报错
shopt -s nullglob
# 获取所有可能的 R1 文件列表
fastq_files=(${RAW_DIR}/*_1.fastq.gz ${RAW_DIR}/*_1.fq.gz)
shopt -u nullglob

count_fastq=${#fastq_files[@]}

if [ "$count_fastq" -eq "0" ]; then
    echo "❌ 错误: 在 ${RAW_DIR} 中未找到 *_1.fastq.gz 或 *_1.fq.gz 文件！"
    exit 1
fi

echo "检测到 $count_fastq 个样本，开始处理..."

for r1_file in "${fastq_files[@]}"
do
    [ -e "$r1_file" ] || continue

    filename=$(basename "${r1_file}")
    
    # --- 关键修改：智能识别后缀 ---
    if [[ "$filename" == *"_1.fastq.gz" ]]; then
        suffix="_1.fastq.gz"
        r2_suffix="_2.fastq.gz"
    elif [[ "$filename" == *"_1.fq.gz" ]]; then
        suffix="_1.fq.gz"
        r2_suffix="_2.fq.gz"
    else
        echo "跳过未知后缀文件: $filename"
        continue
    fi
    
    sample_name=${filename%$suffix}
    r2_file="${RAW_DIR}/${sample_name}${r2_suffix}"
    # -----------------------------

    echo "------------------------------------------"
    echo ">>> 正在处理样本: ${sample_name} (格式: $suffix) <<<"
    echo "------------------------------------------"

    # Step 1: Fastp
    echo "[1/3] 运行 fastp..."
    if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
        if [ ! -f "$r2_file" ]; then
             echo "❌ 错误: 找不到对应的 R2 文件: $r2_file"
             continue
        fi

        fastp -i "${r1_file}" -I "${r2_file}" \
              -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
              -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
              -h "${TRIM_DIR}/${sample_name}_fastp.html" \
              -j "${TRIM_DIR}/${sample_name}_fastp.json" \
              --thread 8 --detect_adapter_for_pe --length_required 25
    else
        echo "   跳过 (文件已存在)"
    fi

    
    # Step 2: Bowtie2 (去rRNA/mtDNA)
    echo "[2/3] 去除 rRNA 和 线粒体..."
    if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
        if [ ! -f "${RRNA_INDEX}.1.bt2" ] && [ ! -f "${RRNA_INDEX}.1.bt2l" ]; then
             echo "错误: rRNA索引未找到: ${RRNA_INDEX}"
             exit 1
        fi

        bowtie2 -p 12 --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" \
                > /dev/null
        
        # 处理 Bowtie2 输出的重命名
        mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
        mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
    else
        echo "   跳过 (文件已存在)"
    fi

    # Step 3: STAR
    echo "[3/3] 运行 STAR..."
    if [ ! -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
        STAR --runThreadN 12 --genomeDir "${STAR_INDEX}" \
             --readFilesIn "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" \
             --readFilesCommand zcat \
             --outFileNamePrefix "${ALIGN_DIR}/${sample_name}." \
             --outSAMtype BAM SortedByCoordinate \
             --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
             --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All \
             --genomeSAsparseD 3

        samtools index "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
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

echo "🎉 所有流程完成。"
echo "所有流程完成。"
