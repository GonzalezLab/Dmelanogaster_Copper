library(chisq.posthoc.test)

# Tolerant
## GIM-012
GIM012 <- read.table("GIM-012.TvsC.tab", row.names = 1, sep = "\t")
GIM012 <- GIM012 %>%
  row_to_names(row_number = 1) 

GIM012<-mutate_all(GIM012, function(x) as.numeric(as.character(x)))

chisq.test(GIM012)
chisq.posthoc.test(GIM012,method = "bonferroni")

## MUN-020
MUN020 <- read.table("MUN-020.TvsC.tab", row.names = 1, sep = "\t")
MUN020 <- MUN020 %>%
  row_to_names(row_number = 1) 

MUN020<-mutate_all(MUN020, function(x) as.numeric(as.character(x)))

chisq.test(MUN020)
chisq.posthoc.test(MUN020,method = "bonferroni")


## MUN-008
MUN008 <- read.table("MUN-008.TvsC.tab", row.names = 1, sep = "\t")
MUN008 <- MUN008 %>%
  row_to_names(row_number = 1) 

MUN008<-mutate_all(MUN008, function(x) as.numeric(as.character(x)))

chisq.test(MUN008)
chisq.posthoc.test(MUN008,method = "bonferroni")


# Sensitive
## AKA-018
AKA018 <- read.table("AKA-018.TvsC.tab", row.names = 1, sep = "\t")
AKA018 <- AKA018 %>%
  row_to_names(row_number = 1) 

AKA018<-mutate_all(AKA018, function(x) as.numeric(as.character(x)))

chisq.test(AKA018)
chisq.posthoc.test(AKA018,method = "bonferroni")

## COR-018
COR018 <- read.table("COR-018.TvsC.tab", row.names = 1, sep = "\t")
COR018 <- COR018 %>%
  row_to_names(row_number = 1) 

COR018<-mutate_all(COR018, function(x) as.numeric(as.character(x)))

chisq.test(COR018)
chisq.posthoc.test(COR018,method = "bonferroni")


## JUT-008
JUT008 <- read.table("JUT-008.TvsC.tab", row.names = 1, sep = "\t")
JUT008 <- JUT008 %>%
  row_to_names(row_number = 1) 

JUT008<-mutate_all(JUT008, function(x) as.numeric(as.character(x)))

chisq.test(JUT008)
chisq.posthoc.test(JUT008,method = "bonferroni")

## TOLERANT: TREATED VS CONTROL
tolerant.TvsC <- read.table("tolerant.TvsC.tab", row.names = 1, sep = "\t")
tolerant.TvsC <- tolerant.TvsC %>%
  row_to_names(row_number = 1) 

tolerant.TvsC<-mutate_all(tolerant.TvsC, function(x) as.numeric(as.character(x)))

chisq.test(tolerant.TvsC)
chisq.posthoc.test(tolerant.TvsC,method = "bonferroni")


## SENSITIVE: TREATED VS CONTROL
sensitive.TvsC <- read.table("sensitive.TvsC.tab", row.names = 1, sep = "\t")
sensitive.TvsC <- sensitive.TvsC %>%
  row_to_names(row_number = 1) 

sensitive.TvsC<-mutate_all(sensitive.TvsC, function(x) as.numeric(as.character(x)))

chisq.test(sensitive.TvsC)
chisq.posthoc.test(sensitive.TvsC,method = "bonferroni")


## TREATED: tolerant VS sensitive
treated.TvsS.tab <- read.table("treated.TvsS.tab", row.names = 1, sep = "\t")
treated.TvsS.tab <- treated.TvsS.tab %>%
  row_to_names(row_number = 1) 

treated.TvsS.tab<-mutate_all(treated.TvsS.tab, function(x) as.numeric(as.character(x)))

chisq.test(treated.TvsS.tab)
chisq.posthoc.test(treated.TvsS.tab,method = "bonferroni")


## CONTROL: tolerant VS sensitive
control.TvsS.tab <- read.table("control.TvsS.tab", row.names = 1, sep = "\t")
control.TvsS.tab <- control.TvsS.tab %>%
  row_to_names(row_number = 1) 

control.TvsS.tab<-mutate_all(control.TvsS.tab, function(x) as.numeric(as.character(x)))

chisq.test(control.TvsS.tab)
chisq.posthoc.test(control.TvsS.tab,method = "bonferroni")


strains_tolerant <- c("GIM-012", "MUN-020", "MUN-008")
strains_sensitive <- c("AKA-018", "COR-018", "JUT-008")
exps <- c("treated","control")

for (exp in exps) {
  for (strain_tolerant in strains_tolerant) {
    for (strain_sensitive in strains_sensitive) {
      contingency_table <- read.table(paste0(strain_tolerant,"_vs_",strain_sensitive,".",exp,".tab"), row.names = 1, sep = "\t")
      contingency_table <- contingency_table %>%
        row_to_names(row_number = 1) 
      
      contingency_table<-mutate_all(contingency_table, function(x) as.numeric(as.character(x)))
      
      print(chisq.test(contingency_table))
      print(chisq.posthoc.test(contingency_table))
      kk<-chisq.posthoc.test(contingency_table)
      kk[kk$Value=="p values",c(3,4)] <- lapply(kk[kk$Value=="p values",c(3,4)], function(x) (p.adjust(x,method="bonferroni",n=54)))
      #write.table(kk,"result_pairs_chisq.posthoc.BF.tab", col.names = T,append = T,quote = F, sep ="\t", row.names = F)
      table_comp<-data.frame(paste0(strain_tolerant," vs ",strain_sensitive," - ",exp),chisq.test(contingency_table)[c(1,3)])
      colnames(table_comp)  <- c("Comparison","X-squared","p-value")
      table_comp$`X-squared`  <- signif(table_comp$`X-squared`, digits=3)
      table_comp$`p-value`  <- signif(table_comp$`p-value`, digits=3)
      table_comp$`p-value` <- p.adjust(table_comp$`p-value`,method="bonferroni",n=9)
      write.table(table_comp,"result_pairs_chisq.BF.tab", col.names = T,append = T,quote = F, sep ="\t", row.names = F)
    }
  }
}
