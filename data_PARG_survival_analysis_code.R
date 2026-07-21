library(survival)
library(ggplot2)
library(survminer)

# ── 1. 读取数据 ───────────────────────────────────────────
data <- read.csv("data_PARG_survival_analysis.csv",
                 stringsAsFactors = FALSE)

cat("样本数：", nrow(data), "\n")
cat("列名：",   paste(colnames(data), collapse = ", "), "\n")

# ── 2. 按 PARG 中位数分组 ────────────────────────────────
median_val  <- median(data$PARG, na.rm = TRUE)
cat("PARG 中位数：", median_val, "\n")

data$Group   <- factor(
  ifelse(data$PARG > median_val, "High", "Low"),
  levels = c("High", "Low")
)
data$OS_time <- data$OS_time / 30  # 天转为月

cat("分组情况：\n")
table(data$Group)

# ── 3. 生存分析 ───────────────────────────────────────────
fit  <- survfit(Surv(OS_time, OS) ~ Group, data = data)
diff <- survdiff(Surv(OS_time, OS) ~ Group, data = data)
p    <- 1 - pchisq(diff$chisq, df = 1)

p_label <- ifelse(p < 0.001, "p < 0.001",
                  paste0("p = ", sprintf("%.3f", p)))

# Cox HR
cox_sum <- summary(coxph(Surv(OS_time, OS) ~ Group, data = data))
HR      <- cox_sum$conf.int[1, "exp(coef)"]
HR_lo   <- cox_sum$conf.int[1, "lower .95"]
HR_hi   <- cox_sum$conf.int[1, "upper .95"]

cat("p 值：", p_label, "\n")
cat("HR =", round(HR, 2),
    "(", round(HR_lo, 2), "-", round(HR_hi, 2), ")\n")
print(fit)

# ── 4. 绘图 ───────────────────────────────────────────────
colors <- c("#F09148", "#427AB2")

surv_plot <- ggsurvplot(
  fit, data,
  conf.int          = TRUE,
  pval              = p_label,
  pval.size         = 5,
  legend.title      = "PARG",
  legend.labs       = c("High", "Low"),
  xlab              = "Time (months)",
  ylab              = "Survival probability",
  break.time.by     = 20,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  palette           = colors,
  surv.median.line  = "hv",
  title             = paste0("PARG Overall Survival (median = ",
                             round(median_val, 3), ")")
)

surv_plot$plot <- surv_plot$plot +
  theme(plot.title = element_text(hjust = 0.5, size = 12,
                                  face = "bold"))

print(surv_plot)

# ── 5. 保存 ───────────────────────────────────────────────
pdf("PARG_survival_median.pdf", width = 8, height = 7)
print(surv_plot)
dev.off()

ggsave("PARG_survival_median.png",
       plot  = surv_plot$plot,
       width = 8, height = 6, dpi = 300)

cat("✓ 已保存\n")
