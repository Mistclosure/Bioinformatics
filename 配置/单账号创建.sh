# 1. 定义变量
NEW_USER="lilin"
INITIAL_PASS="lilin"

# 2. 创建用户 (不交互，直接设置)
# --gecos "" 是为了跳过姓名、电话等详细信息的填写
sudo adduser --disabled-password --gecos "" $NEW_USER

# 3. 设置初始密码为 lilin
echo "$NEW_USER:$INITIAL_PASS" | sudo chpasswd

# 4. 强制用户在下次登录时必须修改密码
sudo chage -d 0 $NEW_USER

# 5. 获取新用户的 UID 并配置资源策略
USER_UID=$(id -u $NEW_USER)
SLICE_NAME="user-${USER_UID}.slice"

echo "正在为新用户 $NEW_USER (UID: $USER_UID) 配置资源策略..."

# 应用 CPU 权重 (100 为基准)
sudo systemctl set-property $SLICE_NAME CPUWeight=100

# 应用内存策略 (最大 250G，忙时保底 60G)
sudo systemctl set-property $SLICE_NAME MemoryMax=250G MemoryLow=60G

# 重新加载配置
sudo systemctl daemon-reload

echo "--------------------------------------------------"
echo "配置完成！"
echo "用户名: $NEW_USER"
echo "初始密码: $INITIAL_PASS (登录后将强制要求修改)"
echo "资源限制: CPUWeight=100, MemoryMax=250G, MemoryLow=60G"
