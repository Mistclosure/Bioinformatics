#!/bin/bash

# ==========================================
# TE Analysis Pipeline (全自动高性能版)
# ==========================================
# 适配配置: 180 Threads / 236GB RAM
# 优化策略: Step0-3 单样本饱和计算 -> Step4 多样本并行计算

# ======================
# 0. 全局资源配置
# ======================
HIGH_THREADS=128  
MID_THREADS=48
LOW_THREADS=16

# 检查环境
for cmd in fastq-dump fasterq-dump fastp bowtie2 STAR TEcount samtools; do
    if ! command -v $cmd &> /dev/null; then
        echo "❌ 错误: 命令 $cmd 未找到！请先运行: conda activate te_env"
        exit 1
    fi
done

# ======================
# 1. 路径配置
# ======================
WORKDIR="/home/qiuzerui/RNA-seq/Phf20_GSE82115"
# 确保目录存在
if [ ! -d "$WORKDIR" ]; then
    echo "❌ 错误: 工作目录 $WORKDIR 不存在！"
    exit 1
fi
cd ${WORKDIR} || exit 1

SRA_DIR="${WORKDIR}/sra"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 注意：请确保这些资源路径是正确的
ANNO_DIR="/home/qiuzerui/RNA-seq/annotations/annotationHv49"
STAR_INDEX="/home/qiuzerui/RNA-seq/indexes/star_index_h38"
GTF_GENE="${ANNO_DIR}/gencode.v49.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/h38_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

mkdir -p ${RAW_DIR} ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}

# ==========================================
# Step 0: SRA 转 FASTQ (内存加速版)
# ==========================================
echo "=== Step 0: 检查 SRA 数据源 ==="

shopt -s nullglob
sra_files=(${SRA_DIR}/*.sra)
shopt -u nullglob

if [ ${#sra_files[@]} -gt 0 ]; then
    echo ">>> 开始 SRA -> FASTQ 转换..."
    for sra_file in "${sra_files[@]}"; do
        filename=$(basename ${sra_file})
        sample_name=${filename%.sra}
        
        # 简单检查是否已存在
        if [ -f "${RAW_DIR}/${sample_name}_1.fastq.gz" ] || [ -f "${RAW_DIR}/${sample_name}_1.fq.gz" ]; then
            continue
        fi

        echo "正在处理: ${sample_name}"
        if command -v fasterq-dump &> /dev/null; then
            fasterq-dump --split-3 -e ${MID_THREADS} -m 2048MB --outdir ${RAW_DIR} --progress ${sra_file}
            
            if command -v pigz &> /dev/null; then
                pigz -p ${MID_THREADS} ${RAW_DIR}/${sample_name}_1.fastq
                pigz -p ${MID_THREADS} ${RAW_DIR}/${sample_name}_2.fastq
            else
                gzip -f ${RAW_DIR}/${sample_name}_1.fastq
                gzip -f ${RAW_DIR}/${sample_name}_2.fastq
            fi
        else
            fastq-dump --split-3 --gzip --outdir ${RAW_DIR} ${sra_file}
        fi
    done
else
    echo "⏭️  SRA 目录无文件或已处理，跳过。"
fi

# ==========================================
# Step 1-3: 预处理和比对循环 (单样本串行，多线程加速)
# ==========================================
echo "=== 开始 Step 1-3: 预处理和比对 ==="

shopt -s nullglob
fastq_files=(${RAW_DIR}/*_1.fastq.gz ${RAW_DIR}/*_1.fq.gz)
shopt -u nullglob

if [ ${#fastq_files[@]} -gt 0 ]; then
    for r1_file in "${fastq_files[@]}"; do
        filename=$(basename "${r1_file}")
        
        # 智能后缀识别
        if [[ "$filename" == *"_1.fastq.gz" ]]; then
            suffix="_1.fastq.gz"
            r2_suffix="_2.fastq.gz"
        elif [[ "$filename" == *"_1.fq.gz" ]]; then
            suffix="_1.fq.gz"
            r2_suffix="_2.fq.gz"
        else 
            continue 
        fi
        
        sample_name=${filename%$suffix}
        r2_file="${RAW_DIR}/${sample_name}${r2_suffix}"

        # 检查是否全部完成，如果bam已存在则彻底跳过
        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            # echo "   -> ${sample_name} 比对已完成，跳过。"
            continue
        fi

        echo ">>> 处理样本: ${sample_name} <<<"

        # Step 1: Fastp
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            echo "[1/3] fastp 质控..."
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25
        fi

        # Step 2: Bowtie2
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            echo "[2/3] Bowtie2 去 rRNA..."
            bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                    -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                    -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                    --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" \
                    > /dev/null 2>&1
            
            mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
            mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
        fi

        # Step 3: STAR
        if [ ! -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "[3/3] STAR 比对..."
            STAR --runThreadN ${HIGH_THREADS} --genomeDir "${STAR_INDEX}" \
                 --readFilesIn "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" \
                 --readFilesCommand zcat \
                 --outFileNamePrefix "${ALIGN_DIR}/${sample_name}." \
                 --outSAMtype BAM SortedByCoordinate \
                 --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
                 --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All \
                 --genomeSAsparseD 3 \
                 --limitBAMsortRAM 60000000000 

            samtools index -@ 32 "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
        fi
    done
else
    echo "⚠️  未检测到原始 FASTQ 文件，假设已进入比对阶段..."
fi

# ==========================================
# Step 4: TEcount (全自动并行定量)
# ==========================================
echo "=== Step 4: TEcount 自动扫描并行定量 ==="

# 1. 启用 nullglob 防止找不到文件报错
shopt -s nullglob
# 自动扫描所有 BAM 文件
bam_files=(${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

count_bams=${#bam_files[@]}

if [ "$count_bams" -eq "0" ]; then
    echo "❌ 错误: 在 ${ALIGN_DIR} 中未找到任何符合要求的 BAM 文件！"
    exit 1
fi

echo "✅ 检测到 $count_bams 个 BAM 文件，准备并行处理..."

# 2. 循环处理每一个检测到的 BAM 文件
for bam_file in "${bam_files[@]}"; do
    
    # --- 关键修改：自动提取样本名 ---
    # basename 命令去掉路径，第二个参数去掉后缀
    sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)

    # 检查输出是否存在
    if [ -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then
        echo "✅ [跳过] ${sample_name} 结果已存在。"
        continue
    fi

    echo "🚀 [后台启动] 正在定量: ${sample_name}"

    # --- 核心命令 (后台运行 &) ---
    (
        TEcount --sortByPos --format BAM --mode multi \
                --GTF "${GTF_GENE}" \
                --TE "${GTF_TE}" \
                --project "${COUNTS_DIR}/${sample_name}" \
                --stranded reverse \
                -b "${bam_file}" \
        && echo "🎉 [完成] ${sample_name}"
    ) & 

done

# 3. 等待所有任务
echo "⏳ 所有任务已投递，正在后台全速计算 (请勿关闭终端)..."
wait
echo "✅ 所有流程圆满结束！"
