# ==========================================
# TE Analysis Pipeline (Modified for Phf20-26.1.23)
# ==========================================
# 硬件: 根据当前系统资源自动调整，默认保留原高配设置
# 请确保在 WSL 或 Linux 环境下运行

# 1. 尝试初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate te_env

# ======================
# 🚀 核心配置区域 (已修改)
# ======================

# 解除 Linux 文件打开数量限制
ulimit -n 65535

# [CPU 策略] (请根据你当前机器的实际核心数适当调整)
HIGH_THREADS=100    
MID_THREADS=80     
LOW_THREADS=50      

# [内存 策略]
DUMP_MEM="8000MB"  # Windows/WSL下可能需要适当降低
STAR_RAM="150000000000" # 约150G，如内存不足请调低

# [路径配置 - 根据你的需求修改]
WORKDIR="/mnt/windowsdata/qiuzerui/Phf20-26.1.23"
RAW_DIR="${WORKDIR}/rawdata"

# 生成的中间文件目录
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

# 注释与索引路径 (根据附件一和描述修改)
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
STAR_INDEX="/mnt/windowsdata/qiuzerui/indexes/star_index_m39"

GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"

# 根据附件一图片，bowtie2索引前缀为 rRNA_mtDNA_index
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# 初始化目录
echo ">>> 正在初始化目录..."
if cd "${WORKDIR}"; then
    echo "✅ 已进入工作目录: $(pwd)"
else
    echo "❌ [报错] 无法进入目录 $WORKDIR，请检查路径是否存在！"
    exit 1
fi
mkdir -p ${TRIM_DIR} ${CLEAN_DIR} ${ALIGN_DIR} ${COUNTS_DIR}

# ==========================================
# Step 0: SRA 转 FASTQ (已跳过)
# ==========================================
# 用户提示 FQ 文件已存在，跳过此步骤
echo "=== Step 0: SRA 转换已跳过 (Raw FASTQ 已就绪) ==="


# ==========================================
# Step 1-3: 智能匹配 & 预处理 & 比对
# ==========================================
echo "=== Step 1-3: 智能匹配模式 (针对 .raw.fastq.gz) ==="

shopt -s nullglob
# 扫描 rawdata 下的所有 gz 文件
all_files=(${RAW_DIR}/*.gz)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    echo "✅ 扫描到 ${#all_files[@]} 个文件，开始处理..."

    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        # --- 智能匹配逻辑 (根据附件二文件名修改) ---
        # 文件名示例: L1MKL2609676-Scr_1_Mixt.R1.raw.fastq.gz
        
        # 1. 如果是 R2 文件，跳过 (我们在处理 R1 时自动找 R2)
        if [[ "$filename" =~ \.R2\.raw\.fastq\.gz$ ]]; then
            continue
        fi

        # 2. 识别 R1 文件并构建 R2 文件名
        if [[ "$filename" =~ \.R1\.raw\.fastq\.gz$ ]]; then
            # 将 .R1.raw.fastq.gz 替换为 .R2.raw.fastq.gz
            r2_filename="${filename/.R1.raw.fastq.gz/.R2.raw.fastq.gz}"
            
            # 提取样本名 (去掉后缀)
            sample_name=$(echo "$filename" | sed 's/\.R1\.raw\.fastq\.gz//')
        else
            # 如果不符合 R1 命名规则，跳过
            continue
        fi

        r2_file="${RAW_DIR}/${r2_filename}"

        # 检查 R2 是否存在
        if [ ! -f "$r2_file" ]; then
            echo "❌ [报错] 样本 $sample_name 缺少 R2 文件 ($r2_filename)，跳过。"
            continue
        fi

        # 检查是否已完成比对 (断点续传)
        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "✅ [跳过] ${sample_name} 比对已完成。"
            continue
        fi

        echo ">>> 正在处理: ${sample_name} <<<"
        echo "    R1: $filename"
        echo "    R2: $r2_filename"

        # [1/3] Fastp
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            echo "   -> [Fastp] 质控..."
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25 2> /dev/null
        fi

        # [2/3] Bowtie2 (去 rRNA)
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            echo "   -> [Bowtie2] 去除 rRNA..."
            # 检查索引是否存在 (检查 .1.bt2 文件)
            if ls "${RRNA_INDEX}"*.bt2* &> /dev/null; then
                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" \
                        > /dev/null 2>&1
                
                # 重命名输出文件 (Bowtie2 输出可能是 .1 或 .1.gz)
                mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
            else
                echo "❌ [报错] rRNA 索引文件未找到 (${RRNA_INDEX}*)，请检查路径！"
                exit 1
            fi
        fi

        # [3/3] STAR (比对)
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
                     --limitBAMsortRAM ${STAR_RAM} > /dev/null

                samtools index -@ ${LOW_THREADS} "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
            else
                 echo "❌ [报错] STAR 索引目录不存在 ($STAR_INDEX)！"
                 exit 1
            fi
        fi
    done
else
    echo "❌ [报错] Rawdata ($RAW_DIR) 目录下未找到任何 .gz 文件！"
fi
# ==========================================
# Step 4: TEcount 定量 (核心逻辑)
# ==========================================
echo "=== Step 4: TEcount 定量 ==="

# 确保在 te_env 环境下运行（如果之前是在此环境下安装的）
# conda activate te_env 

shopt -s nullglob
bam_files=(${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

count_bams=${#bam_files[@]}

if [ "$count_bams" -eq "0" ]; then
    echo "❌ [报错] 未找到 BAM 文件，无法进行定量。"
else
    echo "✅ 准备并行定量 $count_bams 个样本..."
    for bam_file in "${bam_files[@]}"; do
        sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)

        # 断点续传：跳过已完成的样本
        if [ -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then
            echo "✅ [跳过] ${sample_name} 定量已完成。"
            continue
        fi

        echo "🚀 [后台运行] TEcount: ${sample_name}"

        # 使用括号包裹并加 & 符号实现后台并行
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
    
    echo "⏳ 所有任务已投递，正在并行计算中 (请勿关闭终端)..."
    wait
    echo "✅ [Step 4] TEcount 运行结束！"
fi
