# ==========================================
# TE Analysis Pipeline (稳定串行 + featureCounts版)
# ==========================================
# Step 1-3: 串行处理 (节省系统资源)
# Step 4-5: 并行/多线程处理 (加速定量步骤)
# 硬件: AMD EPYC 7R32 (48 Cores) / 250G RAM

# 1. 尝试初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate qiuzerui

# ======================
# 🚀 核心配置区域
# ======================

# --- [关键修复] 解除 Linux 文件打开数量限制 ---
# 默认是 1024，STAR 排序需要更多。提权到 65535。
ulimit -n 65535

# [CPU 策略]
HIGH_THREADS=48   
MID_THREADS=24    
LOW_THREADS=8     

# [内存 策略]
DUMP_MEM="16384MB"
STAR_RAM="100000000000"

# [路径配置]
BASE_DIR="/home/ubuntu/qiuzerui"
WORKDIR="${BASE_DIR}/RNA-seq/Y90C_CMV-Cre"

SRA_DIR="${WORKDIR}/sra"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"

ANNO_DIR="${BASE_DIR}/RNA-seq/annotations/annotationMv38"
STAR_INDEX="${BASE_DIR}/RNA-seq/indexes/star_index_m39"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# 初始化目录
echo ">>> 正在初始化目录..."
if cd "${WORKDIR}"; then
    echo "✅ 已进入: $(pwd)"
else
    echo "❌ [报错] 无法进入目录 $WORKDIR，请检查！"
fi
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
        
        if ls "${RAW_DIR}/${sample_name}"_*.gz &> /dev/null; then
            continue
        fi

        echo "正在处理: ${sample_name}"
        if command -v fasterq-dump &> /dev/null; then
            fasterq-dump --split-3 -e ${MID_THREADS} -m ${DUMP_MEM} --outdir ${RAW_DIR} --progress ${sra_file}
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
echo "=== Step 1-3: 智能匹配模式 (串行处理) ==="

shopt -s nullglob
all_files=(${RAW_DIR}/*.gz)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    echo "✅ 扫描到 ${#all_files[@]} 个文件，开始处理..."

    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        # --- 智能匹配逻辑 ---
        if [[ "$filename" =~ _2\.(fastq|fq)\.gz$ ]] || [[ "$filename" =~ _R2\.(fastq|fq)\.gz$ ]]; then
            continue
        fi

        if [[ "$filename" =~ _R1\.(fastq|fq)\.gz$ ]]; then
            r2_filename="${filename/_R1./_R2.}"
            sample_name=$(echo "$filename" | sed -E 's/_R1\.(fastq|fq)\.gz$//')
        elif [[ "$filename" =~ _1\.(fastq|fq)\.gz$ ]]; then
            r2_filename="${filename/_1./_2.}"
            sample_name=$(echo "$filename" | sed -E 's/_1\.(fastq|fq)\.gz$//')
        else
            continue
        fi

        r2_file="${RAW_DIR}/${r2_filename}"

        if [ ! -f "$r2_file" ]; then
            echo "❌ [报错] 样本 $sample_name 缺少 R2 文件 ($r2_filename)，跳过。"
            continue
        fi

        if [ -f "${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam" ]; then
            echo "✅ [跳过] ${sample_name} 比对已完成。"
            continue
        fi

        echo ">>> 正在处理: ${sample_name} <<<"

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

        # [2/3] Bowtie2
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            echo "   -> [Bowtie2] 去除 rRNA..."
            if ls "${RRNA_INDEX}"* &> /dev/null; then
                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${CLEAN_DIR}/${sample_name}_clean" \
                        > /dev/null 2>&1
                
                mv "${CLEAN_DIR}/${sample_name}_clean.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                mv "${CLEAN_DIR}/${sample_name}_clean.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" 2>/dev/null || mv "${CLEAN_DIR}/${sample_name}_clean.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
            else
                echo "❌ [报错] rRNA 索引未找到，跳过此步！"
            fi
        fi

        # [3/3] STAR
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
    echo "✅ [Step 4] TEcount 运行结束！"
fi

# ==========================================
# Step 5: featureCounts (新增: 多线程全自动)
# ==========================================
echo "=== Step 5: featureCounts 定量 ==="

# 1. 检查软件
if ! command -v featureCounts &> /dev/null; then
    echo "❌ [报错] 未找到 featureCounts！"
    echo "💡 提示: 请运行 conda install -c bioconda subread 安装"
else
    # 2. 定义输出文件 (featureCounts通常输出为一个包含所有样本的大矩阵)
    FC_OUTPUT="${COUNTS_DIR}/all_samples_featureCounts.txt"
    
    # 3. 获取所有 BAM 文件列表
    # (注意：这里使用刚才 Step 4 相同的 bam_files 数组)
    shopt -s nullglob
    bam_files=(${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam)
    shopt -u nullglob
    
    if [ ${#bam_files[@]} -gt 0 ]; then
        if [ -f "${FC_OUTPUT}" ]; then
            echo "✅ [跳过] featureCounts 结果已存在: ${FC_OUTPUT}"
        else
            echo "🚀 [启动] 正在使用 ${HIGH_THREADS} 线程运行 featureCounts..."
            echo "   -> 参数: 双端(-p), 反向链(-s 2), 线程(-T ${HIGH_THREADS})"
            
            # 4. 运行 featureCounts
            # -T: 线程数 (使用 HIGH_THREADS 全速运行)
            # -p: 双端测序
            # -s 2: 反向链特异性 (对应 TEcount --stranded reverse)
            # -a: GTF 注释文件
            # -o: 输出文件路径
            featureCounts -T ${HIGH_THREADS} \
                          -p -s 2 \
                          -a "${GTF_GENE}" \
                          -o "${FC_OUTPUT}" \
                          "${bam_files[@]}" \
                          2>&1 | tee "${FC_OUTPUT}.log"
            
            echo "🎉 [Step 5] featureCounts 完成！结果保存在: ${COUNTS_DIR}"
        fi
    else
        echo "⚠️ 未找到 BAM 文件，跳过 Step 5"
    fi
fi

echo "✅ === 全流程完美结束 ==="
