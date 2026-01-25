#!/bin/bash
set -e

# ==========================================================
# 资源配置参数 (基于 192 线程 / 256G 内存)
# ==========================================================
TARGET_USERS=("zerui" "xiwei" "zlihui")

# 全局用户组上限 (95%): 192 * 0.95 = 182.4 -> 18240%
GLOBAL_CPU_LIMIT="18240%"

# 个人 CPU 上限 (80%): 192 * 0.8 = 153.6 -> 15300%
INDIVIDUAL_CPU_LIMIT="15300%"

# 内存限制
MEM_HIGH="150G"  # 软限制：开始压缩并触发 Swap
MEM_MAX="220G"   # 硬上限：熔断保护，防止宕机

echo "=========================================================="
echo "步骤 1/3：配置系统级 Swap (64G 防护气囊)"
echo "=========================================================="

if [ -f /swapfile ]; then
    echo "✅ /swapfile 已存在，跳过创建。"
else
    echo "正在分配 64G Swap 空间..."
    sudo fallocate -l 64G /swapfile || sudo dd if=/dev/zero of=/swapfile bs=1G count=64
    sudo chmod 600 /swapfile
    sudo mkswap /swapfile
    sudo swapon /swapfile
    echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
    echo "✅ Swap 配置完成。"
fi

echo -e "\n=========================================================="
echo "步骤 2/3：配置全局用户组限制 (确保系统永不假死)"
echo "=========================================================="

# 限制所有普通用户总和不得超过 95% CPU，预留 5% 给系统内核
echo "设置 user.slice 总配额为 $GLOBAL_CPU_LIMIT ..."
sudo systemctl set-property user.slice CPUQuota=$GLOBAL_CPU_LIMIT
echo "✅ 全局防火墙已开启。"

echo -e "\n=========================================================="
echo "步骤 3/3：配置个人资源硬隔离 (CPU 80% + 内存弹性限额)"
echo "=========================================================="

for USERNAME in "${TARGET_USERS[@]}"; do
    if id "$USERNAME" &>/dev/null; then
        USER_UID=$(id -u "$USERNAME")
        SLICE_NAME="user-${USER_UID}.slice"
        
        echo "正在配置用户: $USERNAME (UID: $USER_UID) ..."
        
        # 应用策略：
        # 1. CPUQuota: 单人最高 153 线程，预留 20% 缓冲区
        # 2. MemoryHigh: 150G 开始受压变慢，但不杀进程
        # 3. MemoryMax: 220G 绝对熔断线
        # 4. MemoryLow: 清空旧配置
        sudo systemctl set-property "$SLICE_NAME" \
            CPUQuota=$INDIVIDUAL_CPU_LIMIT \
            MemoryHigh=$MEM_HIGH \
            MemoryMax=$MEM_MAX \
            MemoryLow=
            
        echo "  -> [OK] CPU: $INDIVIDUAL_CPU_LIMIT | 内存: $MEM_HIGH~$MEM_MAX"
    else
        echo "⚠️ 跳过: 找不到用户 '$USERNAME'"
    fi
done

# 强制重载使配置立即生效
sudo systemctl daemon-reload

echo -e "\n=========================================================="
echo "🎉 技术备份完成！服务器已进入【最高防御等级】"
echo "=========================================================="
