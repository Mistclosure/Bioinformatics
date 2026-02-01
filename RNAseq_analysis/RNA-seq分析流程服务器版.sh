# ==========================================
# TE Analysis Pipeline (Phf20-26.1.23 综合版)
# 包含：TEcount(家族) + featureCounts(基因) + TElocal(位点)
# ==========================================

# 1. 尝试初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate te_env

# ======================
# 🚀 核心配置区域 (严格保留原始设置)
# ======================

# 解除 Linux 文件打开数量限制
ulimit -n 65535

# [CPU 策略]
HIGH_THREADS=100    
MID_THREADS=80     
LOW_THREADS=50     

# [内存 策略]
DUMP_MEM="8000MB"  
STAR_RAM="150000000000" 

# [路径配置]
WORKDIR="/mnt/windowsdata/qiuzerui/Phf20-26.1.23"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 注释与索引
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
STAR_INDEX="/mnt/windowsdata/qiuzerui/indexes/star_index_m39"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# --- [TElocal 专用索引] ---
TELOCAL_INDEX="/mnt/windowsdata/qiuzerui/RNAannotations/TElocal/mm39_rmsk_TE.gtf"

# 初始化目录
echo ">>> 正在初始化目录..."
mkdir -p ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}


# ==========================================
# Step 1-3: 智能匹配 & 预处理 & 比对
# ==========================================
echo "=== Step 1-3: 智能匹配模式 (针对 .raw.fastq.gz) ==="

shopt -s nullglob
all_files=(${RAW_DIR}/*.gz)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        if [[ "$filename" =~ \.R2\.raw\.fastq\.gz$ ]]; then continue; fi

        if [[ "$filename" =~ \.R1\.raw\.fastq\.gz$ ]]; then
            r2_filename="${filename/.R1.raw.fastq.gz/.R2.raw.fastq.gz}"
            sample_name=$(echo "$filename" | sed 's/\.R1\.raw\.fastq\.gz//')
        else
            continue
        fi

        r2_file="${RAW_DIR}/${r2_filename}"
        if [ ! -f "$r2_file" ]; then
            echo "❌ [报错] 样本 $sample_name 缺少 R2，跳过。"
            continue
        fi

        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "✅ [跳过] ${sample_name} 比对已完成。"
            continue
        fi

        echo ">>> 正在处理: ${sample_name} <<<"

        # [1/3] Fastp
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25 2> /dev/null
        fi

        # [2/3] Bowtie2
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            if ls "${RRNA_INDEX}"*.bt2* &> /dev/null; then
                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" > /dev/null 2>&1
                
                mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
            fi
        fi

        # [3/3] STAR
        if [ ! -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            STAR --runThreadN ${HIGH_THREADS} --genomeDir "${STAR_INDEX}" \
                 --readFilesIn "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" \
                 --readFilesCommand zcat \
                 --outFileNamePrefix "${ALIGN_DIR}/${sample_name}." \
                 --outSAMtype BAM SortedByCoordinate \
                 --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
                 --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All \
                 --genomeSAsparseD 3 \
                 --limitBAMsortRAM ${STAR_RAM} > /dev/null

            samtools index -@ ${LOW_THREADS} "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
        fi
    done
fi


# ==========================================
# Step 4: TEcount 定量 (保留原始家族水平)
# ==========================================
echo "=== Step 4: TEcount 定量 ==="

shopt -s nullglob
bam_files=(${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

if [ ${#bam_files[@]} -gt 0 ]; then
    for bam_file in "${bam_files[@]}"; do
        sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
        if [ -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then continue; fi

        (
            TEcount --sortByPos --format BAM --mode multi \
                    --GTF "${GTF_GENE}" --TE "${GTF_TE}" \
                    --project "${COUNTS_DIR}/${sample_name}" \
                    --stranded reverse -b "${bam_file}" \
            && echo "🎉 [TEcount 完成] ${sample_name}"
        ) & 
    done
    wait
fi


# ==========================================
# Step 5: featureCounts 定量 (已整合流程)
# ==========================================
echo "=== Step 5: featureCounts 定量 ==="

if ! command -v featureCounts &> /dev/null; then
    echo "❌ [报错] 未找到 featureCounts！"
else
    FC_OUTPUT="${COUNTS_DIR}/all_samples_featureCounts.txt"
    if [ ${#bam_files[@]} -gt 0 ]; then
        if [ -f "${FC_OUTPUT}" ]; then
            echo "✅ [跳过] featureCounts 结果已存在。"
        else
            echo "🚀 [启动] featureCounts (线程: ${HIGH_THREADS})..."
            featureCounts -T ${HIGH_THREADS} \
                          -p -s 2 \
                          -a "${GTF_GENE}" \
                          -o "${FC_OUTPUT}" \
                          "${bam_files[@]}" \
                          2>&1 | tee "${FC_OUTPUT}.log"
            echo "🎉 [Step 5] featureCounts 完成！"
        fi
    fi
fi


# ==========================================
# Step 6: TElocal 定量 (位点水平，新增步骤)
# ==========================================
echo "=== Step 6: TElocal 定量 (Locus Level) ==="

# 内存保护：虽然总线程设为100，但TElocal加载索引极其消耗内存
# 250G 内存环境下，限制 8 个并发样本是兼顾速度与安全的最优解
MAX_JOBS=8 

if [ ${#bam_files[@]} -gt 0 ]; then
    if [ ! -f "${TELOCAL_INDEX}.locInd" ]; then
        echo "❌ [报错] TElocal 索引未找到: ${TELOCAL_INDEX}.locInd"
    else
        for bam_file in "${bam_files[@]}"; do
            sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
            # 使用 _TElocal 后缀防止覆盖 Step 4 的结果
            TELOCAL_PREFIX="${COUNTS_DIR}/${sample_name}_TElocal"

            if [ -f "${TELOCAL_PREFIX}.cntTable" ]; then
                echo "✅ [跳过] ${sample_name} TElocal 已完成。"
                continue
            fi

            (
                TElocal -b "${bam_file}" \
                        --GTF "${GTF_GENE}" \
                        --TE "${TELOCAL_INDEX}" \
                        --project "${TELOCAL_PREFIX}" \
                        --stranded reverse --mode multi \
                && echo "🎉 [TElocal 完成] ${sample_name}"
            ) & 

            # 并行任务管理 (Semaphore)
            while (( $(jobs -r -p | wc -l) >= MAX_JOBS )); do
                sleep 10
            done
        done
        wait
    fi
fi

echo "✅ === 全流程完美结束 (TEcount/featureCounts/TElocal) ==="
