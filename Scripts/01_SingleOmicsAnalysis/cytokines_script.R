#######################################################################
#######               CYTOKINES ANALYSIS MHE DATA               #######
#######################################################################

## By Teresa Rubio
# Raw data are available at supplementary material: S6_Table.xlsx
# Cytokines concentration measurements made by: ELISA
# Cytokines measurements made in: 'Neurological Impairment Laboratory. INCLIVA (Valencia, Spain)'


# Work directory
workDir = "~/Desktop/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines"
setwd(workDir)
options(stringsAsFactors = FALSE)

# Libraries to use
library(readxl)
library(NOISeq)
library(mixOmics)
source(file.path(workDir, "../../auxiliary_functions.R"))


# 1. Data matrix and Metadata ----
  S5_Table = read_excel("~/Desktop/MHERubio2020/MHEarticle/Final_figures/S5_Table.xlsx")
  S5_Table = as.data.frame(S5_Table)
  colnames(S5_Table) = gsub(" \\(../ml\\)", "", colnames(S5_Table))
  
  metadata = S5_Table[, c("SampleName", "Condition")]
  raw_data = t(S5_Table[, -c(1,2)])
  colnames(raw_data) = metadata$SampleName

  #Output
  cytokines_matrix = raw_data
  save(cytokines_matrix, file = file.path(workDir, "cytokines_matrix.Rda"))

# 2. Quality control ----
  # PCA (NOISeq) ----
  myfactors = data.frame(Sample = metadata$SampleName, Type = metadata$Condition)
  exprdata = readData(data=raw_data, factors=myfactors)                   
  myPCA = dat(exprdata, type = "PCA", logtransf = F, norm = T) #Apply log-transformation to improve plots

  plot.scores.gg(myPCA, comp = c(1,2), title = "Score plot - Raw data")
  plot.loadings.gg(myPCA, comp = c(1,2), title = "Loading plot - Raw data")

  # Boxplot ----
  mycolors = c("salmon", "lightblue")[as.factor(metadata$Condition)]
  boxplot(raw_data, las=2, log="y", col=mycolors, names = metadata$SampleName,
          main = "Cytokines raw data", ylab = "log scale")
  legend("topright", y =c("withoutMHE", "whithMHE"), pch = 15, col = c("lightblue", "salmon"))


# 3. Normality and homoscedasticity  ----
  # ANOVA residuals
  mygroup = metadata$Condition
  aov <- apply(raw_data, 1, function(i) {rbind(summary(aov(as.numeric(i) ~ mygroup))[[1]][["F value"]][1],
                                         summary(aov(as.numeric(i) ~ mygroup))[[1]][["Pr(>F)"]][1])})
  rownames(aov) <- c("F value", "Pr(>F)")

  res <- apply(raw_data, 1, function(i) {aov(as.numeric(i) ~ mygroup)$residuals})
  rownames(res) <- colnames(raw_data)

  # qqnorm --> check normality
  par(mar=c(2,2,1,1))
  par(mfrow=c(4,4))
  lapply(1:ncol(res), function(x){qqnorm(res[,x,drop=F], main = colnames(res)[x])})

  # check homocedasticity
  par(mar=c(2,5,1,1))
  par(mfrow=c(4,4))
  lapply(1:nrow(raw_data), function(i){
  plot(x = c(rep(0,5), rep(1,6)), 
       y = c(res[c(1:5),i], res[c(6:11),i]),
       pch = 19,
       cex = 2,
       xlab = "groups", 
       ylab = "residuals",
       xlim = c(-1,2), 
       axes = F,
       main = rownames(raw_data)[i], 
       col = mycolors)
    box()
    axis(2)
    axis(1, at = c(0,1), labels = c("nMHE", "MHE"))
    })  

  
# 4. Wilcoxon test ----
  # wilcox.test() calculates an exact p-value if the samples contain less than 50 finite values and there are no ties.
  # Tied values occur when two or more observations are equal, whether the observations occur in the same sample or in different samples. 
  # In theory, nonparametric tests were developed for continuous distributions where the probability of a tie is zero. In practice, however, ties often occur.
  # In our example, the sample contains less than 50 finite values but there are ties. 
  # That’s why we use wilcox_test in package coin that take into account the presence of ties for calculating p-values.
  
  wilcox_bq = data.frame("Cytokines" = colnames(raw_data), "p.value" = NA, "Z.stat" = NA)
  
  wilcox_bq = apply(raw_data, 1, function(i){
    wilcox = coin::wilcox_test(as.numeric(i) ~ as.factor(mygroup), alternative = "two.sided", mu = 0, paired = FALSE, conf.int = 0.95)
    c("pval" = coin::pvalue(wilcox), "F stat" = coin::statistic(wilcox))
  })
  wilcox_bq = as.data.frame(t(wilcox_bq))
  wilcox_bq$adj.pval = p.adjust(wilcox_bq[,"pval"], method = "BH") #or its alias "fdr"
  
  cytokines_DE = wilcox_bq[wilcox_bq$adj.pval < 0.05,]

  # Plot
  par(mar=c(2,5,1,1))
  par(mfrow=c(4,4))
  meanplot(t(raw_data), mygroup)

  # Output
  save(cytokines_DE, file = file.path(workDir, "cytokines_DE.Rda"))
  save(wilcox_bq, file = "cytokines_DE_all.Rda")

# 5. PLS-DA
  X <- scale(t(raw_data), center = TRUE, scale = TRUE)
  Y <- as.factor(metadata$Condition)             
  
  plsda.res <- plsda(X, Y, ncomp = 8, scale = FALSE) 
  
  # perf --> error rates
  set.seed(2543)
  perf.plsda <- perf(plsda.res, validation = "loo", 
                     progressBar = FALSE, auc = TRUE) 
  
  plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
  
  # plots
  plotIndiv(plsda.res, legend = T, comp = c(1,2), cex = 5)
  plotVar(plsda.res, comp = c(1,2), cex = 5)
  
  plsda.res$loadings$X[order(plsda.res$loadings$X[,1]),][,1,drop=F]

  
# EXTRA: Statistical table to show in the article ----
  median_nMHE <- apply(raw_data[rownames(cytokines_DE), 1:5], 1, median)
  iqr_nMHE <- apply(raw_data[rownames(cytokines_DE), 1:5], 1, IQR)
  
  median_MHE <- apply(raw_data[rownames(cytokines_DE), 6:11], 1, median)
  iqr_MHE <- apply(raw_data[rownames(cytokines_DE), 6:11], 1, IQR)
  
  final <- rbind(median_nMHE, iqr_nMHE, median_MHE, iqr_MHE)
  write.csv(final, file = file.path(workDir, "cytokines_medianIQR.csv"))

