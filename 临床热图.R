
library(ComplexHeatmap)
library(dplyr)
library(grid) 


data <- read.csv("pvc热图.csv")

# 提取需要的临床数据并进行处理
data <- data %>% 
  select(3,5:13)   # 选择需要的列

##排序
data <- data %>%  
  arrange(1,2,3,4,5,6,7,8)

# 定义连续变量 GRB2 的颜色渐变
col_fun_time <- colorRamp2(
  c(3, 1.5, 0),  # 根据值的范围设置
  c("#DC0000FF", "grey", "#1f78b4")  # 从红色到灰色，再到蓝色的渐变
)

ha <- HeatmapAnnotation(
  type = data$group, 
  GRB2 = data$GRB2, 
  level = data$level,
  col = list(
    type = c("tumor" = "#BC3C29FF", "normal" = "#0072B5FF"),
    GRB2 = col_fun_time, # 连续变量的颜色函数
    level = c("0" = "#3C5488FF", "1" = "#FFCCCC", 
              "2" = "#FF6666", "3" = "#CC0000") ))

# 构建zero矩阵
zero_row_data=datarix(nrow=0, ncol=nrow(data))
Hm <- Heatmap(zero_row_data, top_annotation = ha)

draw(Hm, 
     merge_legend = TRUE,                   
     heatmap_legend_side = "bottom",        
     annotation_legend_side = "bottom",     
     width = unit(16, "cm"),                
     height = unit(1, "cm")                 
)

###基因表达热图####
# 读取表达矩阵文件，假设 'gene.txt' 为基因表达数据，行名为基因名
data <- read.csv("二分类热图.csv")

# 按行标准化数据：对每个基因的表达数据进行 log 转换（加 1 防止 0 值），再进行 z-score 标准化
for (i in 1:nrow(data)) {
  data[i, ] <- scale(log(unlist(data[i, ]) + 1, 2))  # 对每一行（基因）进行标准化
}
data <- as.data.frame(data)# 将数据转换为矩阵
data <- na.omit(data)# 去除含有NA的行

# 定义样本分组：假设每组 A、B、C 各有 3 个样本
samples <- rep(c('PVC','Ctrl'), c(3,3))
##排序
data <- data %>%  
  arrange(N1,N2,N3,N4,N5,N6)
#自定义一些颜色
mycolors <- c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd") #藏青-浅绿
mycolors <- c("#073f82", "#1b71b4", "#58a4cf", "#a2cbe3", "#f2f9fe") #藏蓝-浅蓝
mycolors <- c("#eeecdf", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71") #藏蓝-水泥灰
mycolors <- c("#492952", "#82677e", "white", "#59829e", "#1e4668") #脏紫-脏蓝
mycolors <- c("#1f294e", "#5390b5", "#eaebea", "#d56e5e", "#57121d") #经典红-蓝
mycolors <- colorRampPalette(c(c("#6238ff","#ffffff","#ff220e")))(5)

# 绘制热图
heat <- Heatmap(
  data,
  name = "Expression",
  col = mycolors,
  show_column_names =T,#底部显示列名
  cluster_rows = F,  # 关闭基因聚类
  cluster_columns = FALSE,  # 关闭样本聚类
  show_row_names = F  # 不显示基因名称
  
  )

# 展示热图
draw(heat)

# 假设sample_group是一个因子，包含每个样本的分组信息
sample_group <- factor(sample(c('negitive', 'positive'), ncol(data), replace = TRUE))

# 创建列注释
ha_column <- HeatmapAnnotation(
  df = data.frame(Group = sample_group),
  col = list(Group = c(negitive = "blue", positive = "red"))
)

# 创建热图，并添加列注释
heat <- Heatmap(
  data,
  name = "Expression",
  col = mycolors,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)

# 展示热图
draw(heat)

###分类变量
library(ComplexHeatmap)


# 定义离散颜色映射
mycolors <- c("0" = "#339bd3",   # 阴性 → 蓝色
              "1" = "#eb476f")   # 阳性 → 红色

# 绘制热图
heat <- Heatmap(
  data,
  name = "Status",
  col = mycolors,               # ⭐ 关键：离散颜色映射
  show_column_names = TRUE,     # 显示样本名
  cluster_rows = FALSE,         # 不聚类行
  cluster_columns = FALSE,      # 不聚类列
  show_row_names = FALSE,       # 不显示行名（基因名）
  ##column_title = "Binary Status Heatmap (Red=Positive, Blue=Negative)",
  column_title_gp = gpar(fontsize = 10)
)

# 展示热图
draw(heat, heatmap_legend_side = "right")
