library(survival)
library(survminer)
library(data.table)

example <- fread("example_outputKM.tab", header = T, encoding = "UTF-8")

example$strata <- gsub(",","",example$strata)

fit <- survfit(Surv(example$time, example$status) ~ strata, data = example)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = "Dark2")
