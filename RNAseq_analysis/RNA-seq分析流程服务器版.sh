# ==========================================
# TE Analysis Pipeline (Phf20-26.1.23 综合版)
# 包含：TEcount(家族) + featureCounts(基因) + TElocal(位点)
# 修正记录: 修复 TElocal 索引路径报错问题
# ==========================================

# 1. 尝试初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate te_env

# 🚨 预检: 确保 TElocal 已安装
if ! command -v TElocal &> /dev/null; then
    echo "⚠️ [警告] 未找到 TElocal 命令，请确保已运行: pip install TElocal"
fi

# ======================
# 🚀 核心配置区域 (严格保留原始设置)
# ======================

# 解除 Linux 文件打开数量限制
ulimit -n 65535

# [CPU 策略]
HIGH_THREADS=100    
MID_THREADS=80     
LOW_THREADS=50
Featurecount_THREADS=64

# [内存 策略]
DUMP_MEM="8000MB"  
STAR_RAM="150000000000" 

# [路径配置]
WORKDIR="/mnt/disk1/qiuzerui/downloads/Phf8_GSE212779"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 注释与索引
#小鼠注释与索引
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
STAR_INDEX="/mnt/windowsdata/qiuzerui/indexes/star_index_m39"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"
#人类注释与索引
#ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationHv49"
#STAR_INDEX="/mnt/windowsdata/qiuzerui/indexes/star_index_h38"
#GTF_GENE="${ANNO_DIR}/gencode.v49.annotation_PRI.gtf"
#GTF_TE="${ANNO_DIR}/h38_TE.gtf"
#RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# --- [TElocal 专用索引] ---
# ✅ [修正] 必须包含 .locInd 后缀，否则 TElocal 找不到文件
TELOCAL_INDEX="/mnt/windowsdata/qiuzerui/RNAannotations/TElocal/mm39_rmsk_TE.gtf.locInd"

# 初始化目录
echo ">>> 正在初始化目录..."
mkdir -p ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}


# ==========================================
# Step 1-3: 高兼容性匹配模式 (支持 .fq.gz / .fastq.gz)
# ==========================================
echo "=== Step 1-3: 正在扫描 ${RAW_DIR} 下的压缩文件 ==="

shopt -s nullglob
all_files=(${RAW_DIR}/*.gz)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        # 1. 过滤掉所有 R2 标识的文件，确保循环只从 R1 开始
        # 兼容匹配: _2.gz, _R2.gz, .R2.gz, _2.fastq.gz 等
        if [[ "$filename" =~ ([_.]R2[_.]|[_.]2\.) ]]; then
            continue
        fi

        # 2. 核心正则匹配：支持多种 R1 标识和多种后缀名
        # 捕获组说明:
        # ^(.+)            -> ${BASH_REMATCH[1]}: 样本名 (如 SRR21106098)
        # (_1|_R1|\.R1|\.1) -> ${BASH_REMATCH[2]}: R1 标识符
        # (\.f.*)?         -> ${BASH_REMATCH[3]}: 后缀部分 (如 .fastq, .fq, .raw 等，可选)
        # \.gz$            -> 结尾必须是 .gz
        if [[ "$filename" =~ ^(.+)(_1|_R1|\.R1|\.1)(\.f[^.]*)?\.gz$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            r1_id="${BASH_REMATCH[2]}"
            extension="${BASH_REMATCH[3]}"
            
            # 将 R1 标识符中的 '1' 替换为 '2' 来推导 R2 文件名
            r2_id=$(echo "$r1_id" | sed 's/1/2/')
            r2_filename="${sample_name}${r2_id}${extension}.gz"
            r2_file="${RAW_DIR}/${r2_filename}"
        else
            # 如果不是常规的测序数据命名格式，则跳过
            continue
        fi

        # 3. 检查推导出的 R2 文件是否存在
        if [ ! -f "$r2_file" ]; then
            echo "⚠️ [警告] 发现 R1 但未找到对应 R2: $filename (预期 R2: $r2_filename)"
            continue
        fi

        # 4. 检查断点续跑：比对是否已完成
        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "✅ [跳过] ${sample_name} 已完成比对。"
            continue
        fi

        echo ">>> 🚀 正在处理样本: ${sample_name} (格式: ${extension}.gz) <<<"

        # [1/3] Fastp 质控
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25 2> /dev/null
        fi

        # [2/3] Bowtie2 去除 rRNA
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            if ls "${RRNA_INDEX}"*.bt2* &> /dev/null; then
                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" > /dev/null 2>&1
                
                # 处理 bowtie2 输出的带后缀文件名
                mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
            fi
        fi

        # [3/3] STAR 比对
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
            echo "🚀 [启动] featureCounts (线程: ${Featurecount_THREADS})..."
            featureCounts -T ${Featurecount_THREADS} \
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

# 内存保护：限制 8 个并发样本
MAX_JOBS=8 

if [ ${#bam_files[@]} -gt 0 ]; then
    # ✅ [修正] 直接检查变量本身，不要再加 .locInd 后缀
    if [ ! -f "${TELOCAL_INDEX}" ]; then
        echo "❌ [报错] TElocal 索引未找到: ${TELOCAL_INDEX}"
        echo "   请检查文件路径是否正确 (需包含 .locInd)"
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
                        --sortByPos \
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
