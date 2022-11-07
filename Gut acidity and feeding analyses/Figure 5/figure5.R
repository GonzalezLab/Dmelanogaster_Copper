library(reshape2)
library(dplyr)

## fig 5a
dataPH <- read.table("data5A.tab", header = T, sep = "\t")
dataPH$strain <- factor(dataPH$strain,levels=c("GIM-012","MUN-020","MUN-008","AKA-018","COR-018","JUT-008"))
dataPH$x <- factor(dataPH$x, levels=c("Lower pH","Intermediate","Higher pH","No feeding"))
dataPH_treatment <- ggplot(dataPH[dataPH$type=="Treated",], aes(x = x, y = mean.percentage, fill = strain)) + geom_col(color="black",position="dodge") +
  scale_fill_manual(values=color_map) + 
  geom_errorbar(data = dataPH[dataPH$type=="Treated",], aes(x=x, ymin = mean.percentage-sem,
                                                            ymax = mean.percentage+sem,group=strain),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Percentage sampled",title="Gut pH: copper treatment") +coord_cartesian( ylim=c(0,100),xlim=c(0.35,4.65), expand = FALSE ) + theme(
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,size = 18),
    legend.position = "none")
dataPH_control <- ggplot(dataPH[dataPH$type=="Control",], aes(x = x, y = mean.percentage, fill = strain)) + geom_col(color="black",position="dodge") +
  scale_fill_manual(values=color_map) + 
  geom_errorbar(data = dataPH[dataPH$type=="Control",], aes(x=x, ymin = mean.percentage-sem,
                                                            ymax = mean.percentage+sem,group=strain),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Percentage sampled",title="Gut pH: control") +coord_cartesian( ylim=c(0,100),xlim=c(0.35,4.65), expand = FALSE ) + theme(
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,size = 18),
    legend.position = "none")
fig5A<-plot_grid(dataPH_treatment,dataPH_control,nrow = 1)

## fig 5b
dataFeeding <- read.table("data5B.tab", header = T, sep = "\t")
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
    plot.title = element_text(hjust = 0.5,size = 18),
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
    plot.title = element_text(hjust = 0.5,size = 18),
    legend.position = "none")

fig5B<-plot_grid(plot24h,plot40h,nrow = 1)

## fig 5c
dataPH <- read.table("data5C.tab", header = T, sep = "\t")
dataPH$strain <- factor(dataPH$strain,levels=c("w1118","w1118; CG11594"))
color_map2<-c("w1118"="#B69C9B",
              "w1118; CG11594"="#D5DEDD")

dataPH$x <- factor(dataPH$x, levels=c("Lower pH","Intermediate","Higher pH","No feeding"))
dataPH_treatment <- ggplot(dataPH[dataPH$type=="Treated",], aes(x = x, y = mean.percentage, fill = strain)) + geom_col(color="black",position="dodge") +
  scale_fill_manual(values=color_map2) + 
  geom_errorbar(data = dataPH[dataPH$type=="Treated",], aes(x=x, ymin = mean.percentage-sem,
                                                            ymax = mean.percentage+sem,group=strain),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Percentage sampled",title="Gut pH: copper treatment") +coord_cartesian( ylim=c(0,100),xlim=c(0.35,4.65), expand = FALSE ) + theme(
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,size = 18),
    legend.position = "none")
dataPH_control <- ggplot(dataPH[dataPH$type=="Control",], aes(x = x, y = mean.percentage, fill = strain)) + geom_col(color="black",position="dodge") +
  scale_fill_manual(values=color_map2) + 
  geom_errorbar(data = dataPH[dataPH$type=="Control",], aes(x=x, ymin = mean.percentage-sem,
                                                            ymax = mean.percentage+sem,group=strain),
                width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 20) + labs(fill="Strain",x="",y="Percentage sampled",title="Gut pH: control") +coord_cartesian( ylim=c(0,100),xlim=c(0.35,4.65), expand = FALSE ) + theme(
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,size = 18),
    legend.position = "none")
fig5C<-plot_grid(dataPH_treatment,dataPH_control,nrow = 1)

for ( type in unique(dataPH$x) ) {
  for ( exp in c("Control","Treated") ) {
    print(paste(type,exp, t.test(dataPH[dataPH$strain=="w1118; CG11594" & dataPH$type==exp & dataPH$x==type,7:9],dataPH[dataPH$strain=="w1118" & dataPH$type==exp & dataPH$x==type,7:9])$p.value))
  }
}

plot_grid(fig5A,fig5B,fig5C,ncol=1)
ggsave(plot_grid(fig5A,fig5B,fig5C,ncol=1),filename = "fig5.pdf",  device = 'pdf', width = 16, height = 32)
