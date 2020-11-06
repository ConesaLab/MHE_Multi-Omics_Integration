load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_DE.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Metabolomics/metabolomics_DE.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines/cytokines_DE.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines/cytokines_matrix.Rda")

metadata = read.delim("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_patients.txt")
metadata = metadata[,c("SampleName", "Condition")]


# Final table cytokines ----
load("clusters.Rda")
cytMean = t(apply(cytokines_matrix, 1, tapply, metadata$Condition, median))
cytMean = cbind(cytMean, t(apply(cytokines_matrix, 1, tapply, metadata$Condition, IQR)))
cytMean <- data.frame(cytMean [clusters[[1]],])
cytMean$logFC = log2(cytMean$withMHE / cytMean$withoutMHE)
cytMean = round(cytMean, digits = 2)
cyt = data.frame("feature" = rownames(cytMean),
                 "withoutMHE" = paste0(cytMean$withoutMHE, " (", cytMean$withoutMHE.1, ")"),
                 "withMHE" = paste0(cytMean$withMHE, " (", cytMean$withMHE.1, ")"),
                 "logFC" = cytMean$logFC)
write.table(cyt, file = "cyt_logFC.csv", col.names = T, row.names = F, quote = F, sep = "\t", dec = ".")

# Final table genes----
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_matrix.Rda")
chemotaxis = c("CCL5", "CXCL5", "PF4", "PF4V1", "MSMP", "CCR2", "CXCR3", "CMKLR1")
genes <- data.frame(transcriptomics_matrix, check.names = F)[chemotaxis,]
# genes <- 2^genes
genMean = t(apply(genes, 1, tapply, metadata$Condition, median))
genMean = cbind(genMean, t(apply(genes, 1, tapply, metadata$Condition, IQR)))
genMean <- data.frame(genMean)
genMean$logFC = genMean$withMHE - genMean$withoutMHE
genMean = round(genMean, digits = 2)
gen = data.frame("feature" = rownames(genMean),
                 "withoutMHE" = paste0(genMean$withoutMHE, " (", genMean$withoutMHE.1, ")"),
                 "withMHE" = paste0(genMean$withMHE, " (", genMean$withMHE.1, ")"),
                 "logFC" = genMean$logFC)
write.table(gen, file = "gen_logFC.csv", col.names = T, row.names = F, quote = F, sep = "\t", dec = ".")

