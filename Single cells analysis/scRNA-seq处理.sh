#!/bin/bash

# ================= 路径与配置区域 =================
# 1. 原始 SRA 数据来源
SRA_DIR="/mnt/windowsdata/qiuzerui/sra"

# 2. 新的工作与输出目录
OUT_DIR="/mnt/disk1/qiuzerui/coldmouse"

# 3. 参考基因组路径
REF_PATH="/mnt/windowsdata/qiuzerui/scannotations/mouse/refdata-gex-mm10-2020-A"

# 4. 样本列表
SAMPLES=("SRR35688257" "SRR35688258" "SRR35688259" "SRR35688260" "SRR35688261" "SRR35688262")

# ================= 资源配置 (狂暴模式) =================
# 线程数调整为 150，内存 180GB
THREADS=150
MEM_GB=180

# 压缩线程调整为 32 (极速压缩)
ZIP_THREADS=32

# 预先创建新的输出目录并进入
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "当前工作目录已切换至: $(pwd)"
echo "SRA 源文件路径: $SRA_DIR"
echo "配置: 主线程=$THREADS, 压缩线程=$ZIP_THREADS, 内存=${MEM_GB}GB"

# ================= 主程序执行 =================
echo "开始串行处理 6 个样本..."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "--------------------------------------------------------"
    echo ">>>> 正在处理样本: $SAMPLE"
    
    # --- 步骤 1: SRA 转 FASTQ ---
    # 检查目标 FASTQ 是否已存在于新目录下
    if [ ! -f "${SAMPLE}_S1_L001_R1_001.fastq.gz" ]; then
        echo "[$SAMPLE] 步骤 1: 正在将 SRA 转换为 FASTQ (线程: $THREADS)..."
        
        # fasterq-dump 参数说明:
        # -O: 输出到当前目录 ($OUT_DIR)
        # -t: 临时文件也放在当前目录
        fasterq-dump --split-files --include-technical \
                     -e $THREADS \
                     -O . \
                     -t "./tmp_$SAMPLE" \
                     "$SRA_DIR/$SAMPLE"
        
        echo "[$SAMPLE] 步骤 2: 正在进行多线程高速压缩 (每进程 ${ZIP_THREADS} 线程)..."
        # R1 和 R2 并行压缩 (注意：这里总共会占用 2 * ZIP_THREADS 个线程)
        pigz -p $ZIP_THREADS "${SAMPLE}_1.fastq" && mv "${SAMPLE}_1.fastq.gz" "${SAMPLE}_S1_L001_R1_001.fastq.gz" &
        pigz -p $ZIP_THREADS "${SAMPLE}_2.fastq" && mv "${SAMPLE}_2.fastq.gz" "${SAMPLE}_S1_L001_R2_001.fastq.gz" &
        wait
        
        # 清理未压缩的中间文件和临时目录
        rm -f "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq"
        rm -rf "./tmp_$SAMPLE"
    else
        echo "[$SAMPLE] FASTQ 文件已存在于新目录，跳过转换。"
    fi

    # --- 步骤 2: Cell Ranger Count ---
    echo "[$SAMPLE] 步骤 3: 正在运行 Cell Ranger 定量 (线程: $THREADS)..."
    
    if [ ! -d "Output_${SAMPLE}" ]; then
        # 修改点：增加了 --create-bam=false
        cellranger count --id="Output_${SAMPLE}" \
                         --create-bam=false \
                         --transcriptome="$REF_PATH" \
                         --fastqs="$OUT_DIR" \
                         --sample="$SAMPLE" \
                         --localcores=$THREADS \
                         --localmem=$MEM_GB \
                         --jobmode=local
        
        # 【可选】清理原始 FASTQ 以节省空间
        # echo "[$SAMPLE] 清理原始 FASTQ 以节省空间..."
        # rm -f "${SAMPLE}_S1_L001_R1_001.fastq.gz" "${SAMPLE}_S1_L001_R2_001.fastq.gz"
    else
        echo "[$SAMPLE] Output 目录已存在，跳过定量步骤。"
    fi
echo "--------------------------------------------------------"
echo "所有任务完成！结果位于: $OUT_DIR"
