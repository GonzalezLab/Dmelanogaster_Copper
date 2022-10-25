library("AnnotationDbi")
library("org.Dm.eg.db")

mapIds(org.Dm.eg.db,
       keys=c("FBgn0000083","FBgn0000008"), 
       column="SYMBOL", 
       keytype="FLYBASE",
       multiVals="first")  %>% data.frame() %>% rownames_to

mapIds(org.Dm.eg.db,
       keys=c("FBgn0000083","FBgn0000008"), 
       column="ENTREZID", 
       keytype="FLYBASE",
       multiVals="first") %>% as.data.frame()

