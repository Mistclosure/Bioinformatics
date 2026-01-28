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

# ================= 环境变量保险 (防止找不到 cellranger) =================
# 将之前找到的 bin 目录临时加入脚本运行环境
export PATH="/home/zerui/cellranger/bin:$PATH"

# ================= 资源配置 (狂暴模式) =================
# 注意：fasterq-dump 线程过多可能导致 I/O 瓶颈，建议设为 16-32 即可，但这里遵从你的 150
THREADS=150
MEM_GB=180
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
    # 检查目标 FASTQ 是否已存在（检查 R1 即可代表）
    if [ ! -f "${SAMPLE}_S1_L001_R1_001.fastq.gz" ]; then
        echo "[$SAMPLE] 步骤 1: 正在将 SRA 转换为 FASTQ (线程: $THREADS)..."
        
        # fasterq-dump 
        # 警告：如果下载的数据包含 I1 (Index)，这里可能会生成 3 个文件 (_1, _2, _3)
        # 当前脚本假设只有 _1 (R1) 和 _2 (R2)
        fasterq-dump --split-files --include-technical \
                     -e $THREADS \
                     -O . \
                     -t "./tmp_$SAMPLE" \
                     "$SRA_DIR/$SAMPLE"
        
        echo "[$SAMPLE] 步骤 2: 正在进行多线程高速压缩 (每进程 ${ZIP_THREADS} 线程)..."
        
        # R1 和 R2 并行压缩
        pigz -p $ZIP_THREADS "${SAMPLE}_1.fastq" && mv "${SAMPLE}_1.fastq.gz" "${SAMPLE}_S1_L001_R1_001.fastq.gz" &
        pigz -p $ZIP_THREADS "${SAMPLE}_2.fastq" && mv "${SAMPLE}_2.fastq.gz" "${SAMPLE}_S1_L001_R2_001.fastq.gz" &
        wait
        
        # 清理
        rm -f "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq"
        rm -rf "./tmp_$SAMPLE"
    else
        echo "[$SAMPLE] FASTQ 文件已存在于新目录，跳过转换。"
    fi

    # --- 步骤 2: Cell Ranger Count ---
    echo "[$SAMPLE] 步骤 3: 正在运行 Cell Ranger 定量 (线程: $THREADS)..."
    
    if [ ! -d "Output_${SAMPLE}" ]; then
        # 修复点：已加入 --create-bam=false 解决 v8.0 报错
        cellranger count --id="Output_${SAMPLE}" \
                         --create-bam=false \
                         --transcriptome="$REF_PATH" \
                         --fastqs="$OUT_DIR" \
                         --sample="$SAMPLE" \
                         --localcores=$THREADS \
                         --localmem=$MEM_GB \
                         --jobmode=local
    else
        echo "[$SAMPLE] Output 目录已存在，跳过定量步骤。"
    fi
    
    echo "<<<< 样本 $SAMPLE 处理完毕！"

done  # <--- 之前这里缺了 done，一定要加上！

echo "--------------------------------------------------------"
echo "所有任务完成！结果位于: $OUT_DIR"
