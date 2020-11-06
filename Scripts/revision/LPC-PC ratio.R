load("../01_SingleOmicsAnalysis/Metabolomics/metabolomics_matrix.Rda")
load("/home/teresa/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Metabolomics/metabolomics_DE.Rda")

metadata = read.delim("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_patients.txt")
metadata = metadata[,c("SampleName", "Condition")]

met_sig = metabolomics_matrix[rownames(metabolomics_DE),]
pc_sig = met_sig[grep("^PC", rownames(met_sig)),]

ratio_c20.3 = apply(pc_sig, 1, function(x){ x / metabolomics_matrix["lysoPC a C20:3",]})
ratio_c20.4 = apply(pc_sig, 1, function(x){ x / metabolomics_matrix["lysoPC a C20:4",]})


## ANOVA: testar normalidad
group = factor(metadata$Condition)
mytests = apply(t(ratio_c20.4), 1, function (x) aov(x ~ group))
mypvals = sapply(mytests, function (x) anova(x)$`Pr(>F)`[1])
myresiduals <- do.call(rbind, lapply(mytests, function (x) x$residuals))
par(mfrow=c(2,2))
lapply(1:4, function(y){
  qqnorm(myresiduals[y,], main=rownames(myresiduals[y,,drop=F])) 
  qqline(myresiduals[y,])
})
mysd = apply(myresiduals, 1, sd)
mynorm = apply(myresiduals[mysd != 0,], 1, function (x) shapiro.test(x)$p.value)
names(mynorm[mynorm > 0.05])

# t.test
library(limma)
# Build the design matrix for the linear modelling function.
f = factor(metadata$Condition, levels = unique(metadata$Condition))
design = model.matrix(~0 + f)
colnames(design) = c("without_MHE", "with_MHE")

# Create a contrast matrix.
contrast.matrix = makeContrasts("with_MHE-without_MHE", levels = design)

# Apply this contrast matrix to the modeled data and compute statistics for the data.
fit_20.3 = lmFit(t(ratio_c20.3), design)
fit_20.3 = contrasts.fit(fit_20.3, contrast.matrix)
fit_20.3 = eBayes(fit_20.3)

fit_20.4 = lmFit(t(ratio_c20.4), design)
fit_20.4 = contrasts.fit(fit_20.4, contrast.matrix)
fit_20.4 = eBayes(fit_20.4)

ratio_c20.3_DE = topTable(fit_20.3, coef = 1, number = Inf, p.value = "none")
ratio_c20.4_DE = topTable(fit_20.4, coef = 1, number = Inf, p.value = "none")

#output
ratioMean20.3 = data.frame(t(apply(t(ratio_c20.3), 1, tapply, metadata$Condition, mean)))
ratioMean20.4 = data.frame(t(apply(t(ratio_c20.4), 1, tapply, metadata$Condition, mean)))

ratioMean20.3$`limma fdr` = ratio_c20.3_DE[rownames(ratioMean20.3),]$adj.P.Val
ratioMean20.4$`limma fdr` = ratio_c20.4_DE[rownames(ratioMean20.4),]$adj.P.Val


## t.test
p3 = apply(t(ratio_c20.3), 1, function(i){
  test <- t.test(i ~ metadata$Condition)
  test$p.value
})
p4 = apply(t(ratio_c20.4), 1, function(i){
  test <- t.test(i ~ metadata$Condition)
  test$p.value
})
ratioMean20.3$`t-test fdr` <- p.adjust(p = p3, method = "fdr")[rownames(ratioMean20.3)]
ratioMean20.4$`t-test fdr` <- p.adjust(p = p4, method = "fdr")[rownames(ratioMean20.3)]

rownames(ratioMean20.3) = paste(rownames(ratioMean20.3), "/ lysoPC a C20:3")
rownames(ratioMean20.4) = paste(rownames(ratioMean20.4), "/ lysoPC a C20:4")

write.table(round(ratioMean20.3,4), file = "ratio_pc-lyso20_3.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(round(ratioMean20.4,4), file = "ratio_pc-lyso20_4.txt", row.names = T, col.names = T, quote = F, sep = "\t")
