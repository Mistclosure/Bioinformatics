# 定义用户列表
for USERNAME in zerui xiwei lihui shiva; do
    # 获取用户的 UID
    USER_UID=$(id -u $USERNAME)
    
    # 构建 slice 名称 (例如 user-1001.slice)
    SLICE_NAME="user-${USER_UID}.slice"
    
    echo "正在配置用户: $USERNAME (UID: $USER_UID, Slice: $SLICE_NAME)..."
    
    # 应用 CPU 权重 (平分/独占自动调节)
    sudo systemctl set-property $SLICE_NAME CPUWeight=100
    
    # 应用内存策略 (闲时可用满，忙时保底60G)
    # MemoryMax 设为 250G (留一点给系统内核)
    sudo systemctl set-property $SLICE_NAME MemoryMax=250G MemoryLow=60G
    
    echo "配置完成！"
done

# 重新加载配置
sudo systemctl daemon-reload
