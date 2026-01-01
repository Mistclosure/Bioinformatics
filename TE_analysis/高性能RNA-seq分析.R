# ==========================================
# TE Analysis Pipeline (完整通用版)
# ==========================================
# 适配: Ubuntu用户路径 / 智能识别文件名 / 无中断运行

# 1. 尝试初始化 conda (以防直接粘贴时 conda 命令不可用)
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate qiuzerui

# ======================
# 配置区域
# ======================
# 资源
HIGH_THREADS=128  
MID_THREADS=48
LOW_THREADS=16

# 路径 (绝对路径)
BASE_DIR="/home/ubuntu/qiuzerui"
WORKDIR="${BASE_DIR}/RNA-seq/Y90C_CMV-Cre"

# 定义子目录
SRA_DIR="${WORKDIR}/sra"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 参考基因组路径
ANNO_DIR="${BASE_DIR}/RNA-seq/annotations/annotationMv38"
STAR_INDEX="${BASE_DIR}/RNA-seq/indexes/star_index_m39"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# 尝试进入目录
echo ">>> 正在初始化目录..."
if cd "${WORKDIR}"; then
    echo "✅ 已进入: $(pwd)"
else
    echo "❌ [报错] 无法进入目录 $WORKDIR，请检查！"
fi

# 创建文件夹
mkdir -p ${RAW_DIR} ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}


# ==========================================
# Step 0: SRA 转 FASTQ
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
        
        # 检查是否已存在 (支持 .fastq.gz 和 .fq.gz)
        if ls "${RAW_DIR}/${sample_name}"_*.gz &> /dev/null; then
            continue
        fi

        echo "正在处理: ${sample_name}"
        if command -v fasterq-dump &> /dev/null; then
            fasterq-dump --split-3 -e ${MID_THREADS} -m 2048MB --outdir ${RAW_DIR} --progress ${sra_file}
            # 压缩
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
# Step 1-3: 智能匹配 & 预处理 & 比对
# ==========================================
echo "=== Step 1-3: 智能匹配模式 (支持 _1/_R1 和 .fq/.fastq) ==="

shopt -s nullglob
# 扫描所有 .gz 文件
all_files=(${RAW_DIR}/*.gz)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    echo "✅ 扫描到 ${#all_files[@]} 个文件，开始处理..."

    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        # --- 智能匹配逻辑 ---
        # 1. 跳过 R2 文件
        if [[ "$filename" =~ _2\.(fastq|fq)\.gz$ ]] || [[ "$filename" =~ _R2\.(fastq|fq)\.gz$ ]]; then
            continue
        fi

        # 2. 识别 R1 格式
        if [[ "$filename" =~ _R1\.(fastq|fq)\.gz$ ]]; then
            # 格式: Sample_R1.fq.gz
            r2_filename="${filename/_R1./_R2.}"
            sample_name=$(echo "$filename" | sed -E 's/_R1\.(fastq|fq)\.gz$//')
        elif [[ "$filename" =~ _1\.(fastq|fq)\.gz$ ]]; then
            # 格式: Sample_1.fastq.gz
            r2_filename="${filename/_1./_2.}"
            sample_name=$(echo "$filename" | sed -E 's/_1\.(fastq|fq)\.gz$//')
        else
            continue
        fi

        r2_file="${RAW_DIR}/${r2_filename}"

        # 检查 R2 是否存在
        if [ ! -f "$r2_file" ]; then
            echo "❌ [报错] 样本 $sample_name 缺少 R2 文件 ($r2_filename)，跳过。"
            continue
        fi

        # 检查 BAM 是否已存在
        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "✅ [跳过] ${sample_name} 比对已完成。"
            continue
        fi

        echo ">>> 正在处理: ${sample_name} <<<"

        # [1/3] Fastp 质控
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            echo "   -> [Fastp] 质控..."
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25 2> /dev/null
        fi

        # [2/3] Bowtie2 去 rRNA
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            echo "   -> [Bowtie2] 去除 rRNA..."
            if ls "${RRNA_INDEX}"* &> /dev/null; then
                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" \
                        > /dev/null 2>&1
                
                # 重命名输出
                mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
            else
                echo "❌ [报错] rRNA 索引未找到，跳过此步！"
            fi
        fi

        # [3/3] STAR 比对
        if [ ! -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "   -> [STAR] 比对..."
            if [ -d "${STAR_INDEX}" ]; then
                STAR --runThreadN ${HIGH_THREADS} --genomeDir "${STAR_INDEX}" \
                     --readFilesIn "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" \
                     --readFilesCommand zcat \
                     --outFileNamePrefix "${ALIGN_DIR}/${sample_name}." \
                     --outSAMtype BAM SortedByCoordinate \
                     --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
                     --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All \
                     --genomeSAsparseD 3 \
                     --limitBAMsortRAM 60000000000 > /dev/null

                samtools index -@ 32 "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
            else
                 echo "❌ [报错] STAR 索引目录不存在！"
            fi
        fi
    done
else
    echo "❌ [报错] Rawdata 目录下未找到任何 .gz 文件！"
fi


# ==========================================
# Step 4: TEcount (并行定量)
# ==========================================
echo "=== Step 4: TEcount 定量 ==="

shopt -s nullglob
bam_files=(${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

count_bams=${#bam_files[@]}

if [ "$count_bams" -eq "0" ]; then
    echo "❌ [报错] 未找到 BAM 文件，无法进行定量。"
else
    echo "✅ 准备定量 $count_bams 个样本..."
    for bam_file in "${bam_files[@]}"; do
        sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)

        if [ -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then
            echo "✅ [跳过] ${sample_name} 定量已完成。"
            continue
        fi

        echo "🚀 [后台运行] TEcount: ${sample_name}"

        # 后台执行
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
    
    echo "⏳ 所有任务已投递，正在计算中 (请勿关闭终端)..."
    wait
    echo "✅ === 全流程运行结束 ==="
fi
