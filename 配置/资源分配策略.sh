#!/bin/bash

echo "=========================================================="
echo "步骤 1/2：检查并配置 Swap (防爆气囊)"
echo "=========================================================="

# 检查是否已经存在 swapfile
if [ -f /swapfile ]; then
    echo "✅ 检测到 /swapfile 已存在，跳过创建步骤。"
else
    echo "正在创建 64G Swap 文件 (这可能需要几秒钟)..."
    # 创建 64G 空文件
    sudo fallocate -l 64G /swapfile
    # 设置权限 (仅 root 可读写，安全要求)
    sudo chmod 600 /swapfile
    # 格式化为 Swap 分区
    sudo mkswap /swapfile
    # 启用 Swap
    sudo swapon /swapfile
    # 写入开机自动挂载 (防止重启后失效)
    echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
    echo "✅ Swap 创建并启用成功！"
fi

echo -e "\n=========================================================="
echo "步骤 2/2：应用用户资源限制 (CPU硬隔离 + 内存高利用率)"
echo "=========================================================="

# 目标用户列表
TARGET_USERS=("zerui" "xiwei" "zlihui")

for USERNAME in "${TARGET_USERS[@]}"; do
    # 检查用户是否存在，避免报错
    if id "$USERNAME" &>/dev/null; then
        USER_UID=$(id -u "$USERNAME")
        SLICE_NAME="user-${USER_UID}.slice"
        
        echo "正在配置用户: $USERNAME (UID: $USER_UID) ..."
        
        # 核心配置解释：
        # 1. CPUQuota=6000% (60线程): 
        #    硬性锁死。三人加起来 180 线程，永远留 12 线程给系统保命。
        
        # 2. MemoryHigh=150G (软限制/警告线): 
        #    单人可以用到 150G 都是全速。超过后，多余数据进 Swap，速度变慢。
        #    这保证了单人使用时能利用大部分空闲内存。
        
        # 3. MemoryMax=220G (硬上限/熔断线): 
        #    只有极度贪婪（超过220G）才会触发 Kill，防止把服务器物理内存彻底吃光。
        
        # 4. MemoryLow= (清空):
        #    清除旧的保底设置，避免逻辑冲突。
        
        sudo systemctl set-property "$SLICE_NAME" \
            CPUQuota=6000% \
            MemoryHigh=150G \
            MemoryMax=220G \
            MemoryLow=
            
        echo "  -> [OK] CPU锁死60线程 | 内存 150G(软限) - 220G(熔断)"
    else
        echo "⚠️ 跳过: 系统中找不到用户 '$USERNAME'"
    fi
done

# 重载 Systemd 守护进程以确保生效
sudo systemctl daemon-reload

echo -e "\n=========================================================="
echo "🎉 配置全部完成！"
echo "当前策略：单人最大可全速使用 150G 内存，且服务器永远不会因 CPU 争抢而卡死。"
echo "可以使用命令 'systemd-cgtop' 查看实时状态。"
