library("DESeq2")
library("dplyr")
library("tibble")
library("data.table")

###### DE Experiments
# Resistant
# Gene counts
countsFile <- read.table(file="counts_Resistant.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Resistant.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))

colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_resistant_C.vs.T <- res_tb %>%
  dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= log2(1.5))
sigOE_resistant_C.vs.T

sigOE_resistant_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))

# Sensitive
# Gene counts
countsFile <- read.table(file="counts_Sensitive.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Sensitive.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))

colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_sensitive_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_sensitive_C.vs.T

sigOE_sensitive_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))

#########
# Compare Control vs Treated for all strains (no interaction)
# Gene counts
unSortedCounts <- read.table(file="geneCounts.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
countsFile <- unSortedCounts[ , order(names(unSortedCounts))]

# Sample Sheet
unSortedSheet <- read.table(file="sample_sheet.txt", sep="\t", stringsAsFactors=T, header=T)
targetsFile <- unSortedSheet[order(unSortedSheet$Sample),]


rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Strain + Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_all_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_all_C.vs.T

length(intersect(sigOE_all_C.vs.T$gene, c(sigOE_resistant_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene)))
length(intersect(sigOE_all_C.vs.T$gene,sigOE_resistant_C.vs.T$gene))
length(intersect(sigOE_all_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene))
length(intersect(intersect(sigOE_all_C.vs.T$gene,sigOE_resistant_C.vs.T$gene), intersect(sigOE_all_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene)))

#########
# Compare Control vs Treated for all strains (interaction)
# Gene counts
unSortedCounts <- read.table(file="geneCounts.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
countsFile <- unSortedCounts[ , order(names(unSortedCounts))]

# Sample Sheet
unSortedSheet <- read.table(file="sample_sheet.txt", sep="\t", stringsAsFactors=T, header=T)
targetsFile <- unSortedSheet[order(unSortedSheet$Sample),]


rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

colnames(targetsFile) <- c("Sample", "Strain", "Condition", "Batch", "Line")
  
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Line  + Condition + Line:Condition)
ddsFullCountTable



dds <- ddsFullCountTable
dds$Condition <- relevel(dds$Condition, "Control")
dds$Line <- relevel(dds$Line, "Sensitive")

dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds, contrast=c("Condition","Treated","Control"))
res <- results(dds, list( c("Condition_Treated_vs_Control","LineResistant.ConditionTreated") ))
res <- results(dds, list( c("Line_Resistant_vs_Sensitive","LineResistant.ConditionTreated") ))
res <- results (dds, name="LineResistant.ConditionTreated")
res <- results(dds, contrast=c("Line","Resistant","Sensitive"))

head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_all_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_all_C.vs.T
write.table(sigOE_all_C.vs.T, "sigOE_interaction_C.vs.T_interaction.tab", col.names = T, row.names = F, quote = F, sep = "\t")


length(intersect(sigOE_all_C.vs.T$gene, c(sigOE_resistant_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene)))
length(intersect(sigOE_all_C.vs.T$gene,sigOE_resistant_C.vs.T$gene))
length(intersect(sigOE_all_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene))
length(intersect(intersect(sigOE_all_C.vs.T$gene,sigOE_resistant_C.vs.T$gene), intersect(sigOE_all_C.vs.T$gene,sigOE_sensitive_C.vs.T$gene)))

write.table(sigOE_all_C.vs.T, "sigOE_all_C.vs.T.tab", col.names = T, row.names = F, quote = F, sep = "\t")


#########
# Compare Resistant vs Sensitive in control
# Gene counts
countsFile <- read.table(file="counts_Controls.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")

# Sample Sheet
targetsFile <- read.table(file="sample_sheet_Controls.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample


ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Resistance)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Resistance <- relevel(dds$Resistance, "Sensitive")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_control_S.vs.R <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_control_S.vs.R

length(intersect(sigOE_all_S.vs.R$gene, sigOE_control_S.vs.R$gene))

#######
# Treated table
geneCounts <- read.table("geneCounts.txt", header = T)
treated <- fread("sample_sheet_Treated.txt", header = T)
geneCountsTreated <- geneCounts[, c("Geneid", treated$Sample)]
write.table(geneCountsTreated, "counts_Treated.txt", col.names = T, sep = "\t", quote = F, row.names = F)

# Compare Resistant vs Sensitive in treatment
# Gene counts
unSortedCounts <- read.table(file="counts_Treated.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
countsFile <- unSortedCounts[ , order(names(unSortedCounts))]

# Sample Sheet
unSortedSheet <- read.table(file="sample_sheet_Treated.txt", sep="\t", stringsAsFactors=T, header=T)
targetsFile <- unSortedSheet[order(unSortedSheet$Sample),]

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Resistance)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Resistance <- relevel(dds$Resistance, "Sensitive")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_treated_S.vs.R <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_treated_S.vs.R

length(intersect(sigOE_control_S.vs.R$gene, sigOE_treated_S.vs.R$gene))


# ANALYSIS PER STRAIN

#MUN33
# Gene counts
countsFile <- read.table(file="counts_M3_MUN33.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_MUN33.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_MUN33_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_MUN33_C.vs.T

sigOE_MUN33_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
         DownRegulat = sum(log2FoldChange<0))

length(intersect(sigOE_control_S.vs.R$gene, sigOE_treated_S.vs.R$gene))

sigOE_MUN33_C.vs.T  %>%  
  mutate(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))
sigOE_MUN33_C.vs.T$Regulation <- ifelse(sigOE_MUN33_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_MUN33_C.vs.T[,c("gene","Regulation")], "sigOE_MUN33_C.vs.T.lst", quote = F, row.names = F, col.names = F)

#FIN19
# Gene counts
countsFile <- read.table(file="counts_F9_FIN19.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_FIN19.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_FIN19_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_FIN19_C.vs.T

sigOE_FIN19_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))
sigOE_FIN19_C.vs.T$Regulation <- ifelse(sigOE_FIN19_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_FIN19_C.vs.T[,c("gene","Regulation")], "sigOE_FIN19_C.vs.T.lst", quote = F, row.names = F, col.names = F)

#BZ37
# Gene counts
countsFile <- read.table(file="counts_B7_BZ37.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_BZ37.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_BZ37_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_BZ37_C.vs.T

sigOE_BZ37_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))

sigOE_BZ37_C.vs.T$Regulation <- ifelse(sigOE_BZ37_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_BZ37_C.vs.T[,c("gene","Regulation")], "sigOE_BZ37_C.vs.T.lst", quote = F, row.names = F, col.names = F)

#DEN62
# Gene counts
countsFile <- read.table(file="counts_D6_DEN62.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_DEN62.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)

res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_DEN62_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_DEN62_C.vs.T

sigOE_DEN62_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))

sigOE_DEN62_C.vs.T$Regulation <- ifelse(sigOE_DEN62_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_DEN62_C.vs.T[,c("gene","Regulation")], "sigOE_DEN62_C.vs.T.lst", quote = F, row.names = F, col.names = F)

#GIM44
# Gene counts
countsFile <- read.table(file="counts_G4_GIM44.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_GIM44.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_GIM44_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_GIM44_C.vs.T

sigOE_GIM44_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))
sigOE_GIM44_C.vs.T$Regulation <- ifelse(sigOE_GIM44_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_GIM44_C.vs.T[,c("gene","Regulation")], "sigOE_GIM44_C.vs.T.lst", quote = F, row.names = F, col.names = F)

#MUN24
# Gene counts
countsFile <- read.table(file="counts_M4_MUN24.txt", sep="\t", stringsAsFactors=T, header=T, row.names="Geneid")
# Sample Sheet
targetsFile <- read.table(file="sample_sheet_MUN24.txt", sep="\t", stringsAsFactors=T, header=T)

rownames(targetsFile) <- targetsFile$Sample
identical(colnames(countsFile), rownames(targetsFile))
colnames(countsFile) <- targetsFile$Sample

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsFile,
  colData = targetsFile,
  design = ~ Treatment)
ddsFullCountTable

dds <- ddsFullCountTable
dds$Treatment <- relevel(dds$Treatment, "Control")
dds <- DESeq(dds)
plotDispEsts(dds)

res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
head(resOrdered)
resultsNames(dds)


res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE_MUN24_C.vs.T <- res_tb %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
sigOE_MUN24_C.vs.T

sigOE_MUN24_C.vs.T  %>%  
  summarize(UpRegulated = sum(log2FoldChange>0),
            DownRegulat = sum(log2FoldChange<0))

sigOE_MUN24_C.vs.T$Regulation <- ifelse(sigOE_MUN24_C.vs.T$log2FoldChange > 0, "UPREGULATED", "DOWNREGULATED")
write.table(sigOE_MUN24_C.vs.T[,c("gene","Regulation")], "sigOE_MUN24_C.vs.T.lst", quote = F, row.names = F, col.names = F)


### PLOTS
## clusterProfiler
library(clusterProfiler)
organism = "org.Dm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

allGenes <- res_tb$gene

allGenes_symbol<- select(org.Dm.eg.db, 
                         keys = allGenes,
                         columns = c("SYMBOL"),
                         keytype = "FLYBASE")

sigOE_all_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                         keys = sigOE_all_C.vs.T$gene,
                                         columns = c("SYMBOL"),
                                         keytype = "FLYBASE")

allGenes_symbol<- select(org.Dm.eg.db, 
                           keys = allGenes,
                           columns = c("SYMBOL"),
                           keytype = "FLYBASE")

result_all_C.vs.T <- enrichGO(sigOE_all_C.vs.T_genes_symbol$SYMBOL, 
                              universe=allGenes_symbol$SYMBOL,
                              OrgDb=organism,
                              keyType ="SYMBOL",
                              ont="ALL",
                              qvalueCutoff=0.01,
                              pvalueCutoff=0.05)

ggsave("../GO/result_all_C.vs.T.png",
       dotplot(result_all_C.vs.T, showCategory=15, title="All lines (control vs. treatment)"),
       width=14, height=8)
write.table(result_all_C.vs.T, file = "../GO/clusterProfiler_table_result_all_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)


sigOE_control_S.vs.R_genes_symbol<- select(org.Dm.eg.db, 
                                       keys = sigOE_control_S.vs.R$gene,
                                       columns = c("SYMBOL"),
                                       keytype = "FLYBASE")

result_control_S.vs.R <- enrichGO(sigOE_control_S.vs.R_genes_symbol$SYMBOL, 
                                  universe=allGenes_symbol$SYMBOL,
                                  OrgDb=organism,
                                  keyType ="SYMBOL",
                                  ont="ALL",
                                  qvalueCutoff=0.01,
                                  pvalueCutoff=0.05)

ggsave("../GO/result_control_S.vs.R.png",
       dotplot(result_control_S.vs.R, showCategory=15, title="Control (sensitive vs. tolerant)"),
       width=14, height=8)

write.table(result_control_S.vs.R, file = "../GO/clusterProfiler_table_result_control_S.vs.R.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)


sigOE_treated_S.vs.R_genes_symbol<- select(org.Dm.eg.db, 
                                           keys = sigOE_treated_S.vs.R$gene,
                                           columns = c("SYMBOL"),
                                           keytype = "FLYBASE")

result_treated_S.vs.R <- enrichGO(sigOE_treated_S.vs.R_genes_symbol$SYMBOL, 
                                  universe=allGenes_symbol$SYMBOL,
                                  OrgDb=organism,
                                  keyType ="SYMBOL",
                                  ont="ALL",
                                  qvalueCutoff=0.01,
                                  pvalueCutoff=0.05)

ggsave("../GO/result_treated_S.vs.R.png",
       dotplot(result_treated_S.vs.R, showCategory=15, title="Treatment (sensitive vs. tolerant)"),
       width=14, height=8)

write.table(result_treated_S.vs.R, file = "../GO/clusterProfiler_table_result_treated_S.vs.R.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)


sigOE_DEN62_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                       keys = sigOE_DEN62_C.vs.T$gene,
                                       columns = c("SYMBOL"),
                                       keytype = "FLYBASE")

result_DEN62_C.vs.T <- enrichGO(sigOE_DEN62_C.vs.T_genes_symbol$SYMBOL, 
                              universe=allGenes_symbol$SYMBOL,
                              OrgDb=organism,
                              keyType ="SYMBOL",
                              ont="ALL",
                              qvalueCutoff=0.01,
                              pvalueCutoff=0.05)

ggsave("../GO/result_DEN62_C.vs.T.png",
       dotplot(result_DEN62_C.vs.T, showCategory=15, title="JUT-008 (control vs. treatment)"),
       width=14, height=8)
write.table(result_DEN62_C.vs.T, file = "../GO/clusterProfiler_table_result_DEN62_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)



sigOE_BZ37_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                         keys = sigOE_BZ37_C.vs.T$gene,
                                         columns = c("SYMBOL"),
                                         keytype = "FLYBASE")

result_BZ37_C.vs.T <- enrichGO(sigOE_BZ37_C.vs.T_genes_symbol$SYMBOL, 
                                universe=allGenes_symbol$SYMBOL,
                                OrgDb=organism,
                                keyType ="SYMBOL",
                                ont="ALL",
                                qvalueCutoff=0.01,
                                pvalueCutoff=0.05)

ggsave("../GO/result_BZ37_C.vs.T.png",
       dotplot(result_BZ37_C.vs.T, showCategory=15, title="COR-018 (control vs. treatment)"),
       width=14, height=8)
write.table(result_BZ37_C.vs.T, file = "../GO/clusterProfiler_table_result_BZ37_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)



sigOE_FIN19_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                        keys = sigOE_FIN19_C.vs.T$gene,
                                        columns = c("SYMBOL"),
                                        keytype = "FLYBASE")

result_FIN19_C.vs.T <- enrichGO(sigOE_FIN19_C.vs.T_genes_symbol$SYMBOL, 
                               universe=allGenes_symbol$SYMBOL,
                               OrgDb=organism,
                               keyType ="SYMBOL",
                               ont="ALL",
                               qvalueCutoff=0.01,
                               pvalueCutoff=0.05)

ggsave("../GO/result_FIN19_C.vs.T.png",
       dotplot(result_FIN19_C.vs.T, showCategory=15, title="AKA-008 (control vs. treatment)"),
       width=14, height=8)
write.table(result_FIN19_C.vs.T, file = "../GO/clusterProfiler_table_result_FIN19_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)


sigOE_GIM44_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                         keys = sigOE_GIM44_C.vs.T$gene,
                                         columns = c("SYMBOL"),
                                         keytype = "FLYBASE")

result_GIM44_C.vs.T <- enrichGO(sigOE_GIM44_C.vs.T_genes_symbol$SYMBOL, 
                                universe=allGenes_symbol$SYMBOL,
                                OrgDb=organism,
                                keyType ="SYMBOL",
                                ont="ALL",
                                qvalueCutoff=0.01,
                                pvalueCutoff=0.05)

ggsave("../GO/result_GIM44_C.vs.T.png",
       dotplot(result_GIM44_C.vs.T, showCategory=15, title="GIM-012 (control vs. treatment)"),
       width=14, height=8)
write.table(result_GIM44_C.vs.T, file = "../GO/clusterProfiler_table_result_GIM44_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)



sigOE_MUN24_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                         keys = sigOE_MUN24_C.vs.T$gene,
                                         columns = c("SYMBOL"),
                                         keytype = "FLYBASE")

result_MUN24_C.vs.T <- enrichGO(sigOE_MUN24_C.vs.T_genes_symbol$SYMBOL, 
                                universe=allGenes_symbol$SYMBOL,
                                OrgDb=organism,
                                keyType ="SYMBOL",
                                ont="ALL",
                                qvalueCutoff=0.01,
                                pvalueCutoff=0.05)

ggsave("../GO/result_MUN24_C.vs.T.png",
       dotplot(result_MUN24_C.vs.T, showCategory=15, title="MUN-008 (control vs. treatment)"),
       width=14, height=8)
write.table(result_MUN24_C.vs.T, file = "../GO/clusterProfiler_table_result_MUN24_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)





sigOE_MUN33_C.vs.T_genes_symbol<- select(org.Dm.eg.db, 
                                         keys = sigOE_MUN33_C.vs.T$gene,
                                         columns = c("SYMBOL"),
                                         keytype = "FLYBASE")

result_MUN33_C.vs.T <- enrichGO(sigOE_MUN33_C.vs.T_genes_symbol$SYMBOL, 
                                universe=allGenes_symbol$SYMBOL,
                                OrgDb=organism,
                                keyType ="SYMBOL",
                                ont="ALL",
                                qvalueCutoff=0.01,
                                pvalueCutoff=0.05)

ggsave("../GO/result_MUN33_C.vs.T.png",
       dotplot(result_MUN33_C.vs.T, showCategory=15, title="MUN-020 (control vs. treatment)"),
       width=14, height=8)
write.table(result_MUN33_C.vs.T, file = "../GO/clusterProfiler_table_result_MUN33_C.vs.T.tab", quote = F, sep = "\t", row.names = F,
            col.names = T)


plot_grid(dotplot(result_GIM44_C.vs.T, showCategory=15, title="GIM-012 (control vs. treatment)"),
          dotplot(result_MUN24_C.vs.T, showCategory=15, title="MUN-008 (control vs. treatment)"),
          dotplot(result_MUN33_C.vs.T, showCategory=15, title="MUN-020 (control vs. treatment)"),
          dotplot(result_BZ37_C.vs.T, showCategory=15, title="COR-018 (control vs. treatment)"),
          dotplot(result_DEN62_C.vs.T, showCategory=15, title="JUT-008 (control vs. treatment)"),
          dotplot(result_FIN19_C.vs.T, showCategory=15, title="AKA-008 (control vs. treatment)"),
          labels = "AUTO", label_size = 12, ncol=3)

ggsave("../GO/result_individualStrains_C.vs.T.png",
       plot_grid(dotplot(result_GIM44_C.vs.T,  showCategory=15, title="GIM-012 (control vs. treatment)"),
                 dotplot(result_MUN24_C.vs.T, showCategory=15, title="MUN-008 (control vs. treatment)"),
                 dotplot(result_MUN33_C.vs.T, showCategory=15, title="MUN-020 (control vs. treatment)"),
                 dotplot(result_BZ37_C.vs.T, showCategory=15, title="COR-018 (control vs. treatment)"),
                 dotplot(result_DEN62_C.vs.T, showCategory=15, title="JUT-008 (control vs. treatment)"),
                 dotplot(result_FIN19_C.vs.T, showCategory=15, title="AKA-008 (control vs. treatment)"),
                 labels = "AUTO", label_size = 24, ncol=3, align = "hv"), width=28, height=12)


ggsave("../GO/result_all_S.vs.R.png",
       plot_grid(dotplot(result_control_S.vs.R,  showCategory=15, title="Control (sensitive vs. tolerant)"),
                 dotplot(result_treated_S.vs.R, showCategory=15, title="Treatment (sensitive vs. tolerant)"),
                 labels = "AUTO", label_size = 24, ncol=2, align = "hv"), width=16, height=8)

