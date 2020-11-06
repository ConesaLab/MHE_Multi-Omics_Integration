load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_DE.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_matrix.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines/cytokines_DE_all.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines/cytokines_matrix.Rda")

metadata = read.delim("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_patients.txt")
metadata = metadata[,c("SampleName", "Condition")]


genes <- c("IL6",   "IL18",  "IL4",   "IL10",  "IL12A",  "IL13",  "IL15",
           "IL17A",  "IL22",  "CXCL13", "CCL20",  "CX3CL1", "TNF",  "TGFB1")

# cytokines info
cytMean = t(apply(cytokines_matrix, 1, tapply, metadata$Condition, median))
fc = cytMean[,1] / cytMean[,2]
logFC = log2(fc)

table_cyt <- data.frame("cytFDR" = wilcox_bq$adj.pval, "cytLOGFC" = logFC)
table_cyt$name = genes

# messenger info
mygenes = transcriptomics_matrix[rownames(transcriptomics_matrix) %in% genes,]
logfc2 = transcriptomics_DE[genes,"logFC",drop=F]

table_gene <- transcriptomics_DE[rownames(mygenes),c("adj.P.Val", "logFC")]
colnames(table_gene) = c("genFDR", "genLOGFC")

table_gene$name = rownames(table_gene)


# final table
table <- merge(table_cyt, table_gene, by.x = "name", by.y ="name", all = T)

b = round(table[-1], digits = 2)
table <- cbind(table$name, b)

xlsx::write.xlsx(table, file = "Table_S4_final.xlsx")



