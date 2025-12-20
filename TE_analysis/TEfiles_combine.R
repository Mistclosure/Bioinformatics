# 1. 获取工作路径
#请修改变量名称
wd <- '//wsl.localhost/Ubuntu/home/qiuzerui/BRPF1/'
filesname <- '_counts.csv'
#############
setwd(wd)
library(data.table)

# 获取所有 .cntTable 文件
files <- list.files(paste0(wd,'/counts'), 
                    pattern="\\.cntTable$", recursive=TRUE, full.names=TRUE)

# 2. 读取并整合
te_list <- lapply(files, function(f){
  dt <- fread(f)
  sample <- gsub(".*/(.*)\\.cntTable$", "\\1", f)
  colnames(dt) <- c("RepeatID", sample)
  return(dt)
})

# 3. 合并为一个大表
te_merged <- Reduce(function(x, y) merge(x, y, by="RepeatID", all=TRUE), te_list)
te_merged[is.na(te_merged)] <- 0  # 将NA填0

# 4. 输出文件
write.csv(te_merged, filesname, quote=FALSE,row.names = F)
cat("✅ 已生成 counts/ERV_counts.TE_count.txt\n")
