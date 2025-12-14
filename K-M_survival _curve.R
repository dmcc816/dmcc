
####K-M生存曲线####
#读取生存数据，数据框包含OS，OS_time，gene

library(survival)
library("ggplot2")
library("survminer")

data <- read.csv("蛋白组生存曲线.csv",sep = ",",header = T)
col_need<- c("Group","OS","OS_time")
data <- data[,col_need]

##计算最佳cutpoint
cutpoint <- surv_cutpoint(data = data,
                         time = "OS_time",  
                         event = "OS",  
                         variables = c("PSMC2"))  # 需要计算的数据列名，可多个

# 提取 X 的最佳 cutpoint
cutpoint <- cutpoint$cutpoint["PSMC2", "cutpoint"]

# 用 cutpoint 分组（替代 median）
data$Group <- ifelse(data$PSMC2 > cutpoint, "high", "low")

write.csv(data,"data.csv")

data$OS_time <- data$OS_time / 30#month,按照数据修改

fit <- survfit(Surv(OS_time, OS)~Group,data)#risk可以自定义
diff=survdiff(Surv(OS_time, OS)~Group,data = data)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

# 自定义颜色向量
colors <- c("#CC99FF","#0099CC")  # 分别对应 "high" 和 "low","red","yellow"

#绘制训练集的K-M曲线
surv_plot <- ggsurvplot(fit, data,
                        conf.int = T,
                        pval = pValue,
                        pval.size = 6,
                        legend.title = "PSMC2",
                        legend.labs = c("high","low"),##,"-+","--"
                        #xlab = "Time (months)",
                        break.time.by = 20,
                        risk.table = F,
                        risk.table.height = 0.3,
                        palette = colors ,
                        surv.median.line = "hv"
)

# 修改生存曲线图的标题
surv_plot$plot <- surv_plot$plot + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust =0, size = 16, face = "bold"))

# 打印最终的图表
print(surv_plot)


