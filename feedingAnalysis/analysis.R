library(reshape2)
library(dplyr)
dataFeeding <- read.table("data.tab", header = T, sep = "\t")
dataFeeding$strain_no<-NULL
dataFeeding <- melt(dataFeeding, id.vars=c("strain","time","status"))

dataFeedingFCsample <- dataFeeding %>% 
  group_by(strain,time,variable) %>% 
  dplyr::summarize(
    FC = value[status == "treated"]/value[status == "control"]
  ) %>% 
  ungroup()


tmp<-dataFeeding %>% group_by(strain, time, status) %>% dplyr::summarize(mean=mean(value))
FC<-tmp %>% 
  group_by(strain,time) %>% 
  dplyr::summarize(
    FC = mean[status == "treated"]/mean[status == "control"]
  ) %>% 
  ungroup()
sd_controls <- dataFeeding %>% filter(status=="control") %>% group_by(strain,time) %>% dplyr::summarize(stdev=sd(value)/mean(value))
FC<-merge(FC,sd_controls,by=c("strain","time"))

color_map <- c("GIM-012"="#E62D23",
               "MUN-020"="#F18737",
               "MUN-008"="#FDD241",
               "AKA-018"="#C8D441",
               "COR-018"="#B5D9D3",
               "JUT-008"="#72BDD5")
FC$strain <- factor(FC$strain,levels=c("GIM-012","MUN-020","MUN-008","AKA-018","COR-018","JUT-008"))


plot24h<- ggplot(FC[FC$time=="24h",], aes(x = strain, y = FC, fill = strain)) + geom_col(color="black") +
  scale_fill_manual(values=color_map) + 
  geom_errorbar(data = FC[FC$time=="24h",], aes(x=strain, ymin = FC-(FC*stdev),
                               ymax = FC+(FC*stdev)),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Fold-change consumption",title="Feeding: 24h copper treatment") +coord_cartesian( ylim=c(0,1.1),xlim=c(-0.1,7), expand = FALSE ) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none")


plot40h<-ggplot(FC[FC$time=="40h",], aes(x = strain, y = FC, fill = strain)) + geom_col(color="black") +
  scale_fill_manual(values=color_map) + 
  geom_errorbar(data = FC[FC$time=="40h",], aes(x=strain, ymin = FC-(FC*stdev),
                                                ymax = FC+(FC*stdev)),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Fold-change consumption",title="Feeding: 40h copper treatment") +coord_cartesian( ylim=c(0,1.1),xlim=c(-0.1,7), expand = FALSE ) + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none")

plot_grid(plot24h,plot40h,nrow = 1)
ggsave(plot_grid(plot24h,plot40h,nrow = 1),filename = "fig5B.pdf",  device = 'pdf', width = 16)

strains_tolerant <- c("GIM-012", "MUN-020", "MUN-008")
strains_sensitive <- c("AKA-018", "COR-018", "JUT-008")
exps <- c("treated","control")

exps <- c("40h","24h")
strains_tolerant <- c("GIM-012", "MUN-020", "MUN-008")
strains_sensitive <- c("AKA-018", "COR-018", "JUT-008")
for (exp in exps) {
  for (strain_tolerant in strains_tolerant) {
    for (strain_sensitive in strains_sensitive) {
      pvalue_corrected <- p.adjust(t.test(dataFeedingFCsample[dataFeedingFCsample$time==exp & dataFeedingFCsample$strain==strain_tolerant,]$FC,
                                          dataFeedingFCsample[dataFeedingFCsample$time==exp & dataFeedingFCsample$strain==strain_sensitive,]$FC)$p.value,method = "bonferroni",n=9)
      write.table(paste(strain_tolerant,strain_sensitive,exp,pvalue_corrected,sep="\t"),file="result_comparison.tab",col.names = F, row.names = F, quote=F,sep="\t",append=T)
    }
  }
}
      
