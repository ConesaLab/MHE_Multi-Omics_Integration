#######################################################################
#######             MULTIOMICS ANALYSIS - MHE DATA              #######
#######################################################################

## By Teresa Rubio
# Input data is the output of MHERubio2020/Scripts/01_SingleOmicsAnalysis:
#     Transcriptomics: transcriptomics_matrix.Rda & transcriptomics_DE.Rda
#     Metabolomics: metabolomics_matrix.Rda & metabolomics_DE.Rda
#     Cytokines: cytokines_matrix.Rda & cytokines_DE.Rda
# 3 main steps in this analysis: 
#     1) Clustering: metabolomics + cytokines.
#     2) PLS performance: multiple PLS models.
#     3) Enrichment analysis: biological interpretation of each PLS model.


# Work directory and metadata
workDir = "~/Desktop/MHERubio2020/Scripts/02_MultiOmicsAnalysis"
setwd(workDir)
load(file.path(workDir, "metadata.Rda"), verbose = TRUE)
options(stringsAsFactors = FALSE)

# Input directory and data
dataDir = "~/Desktop/MHERubio2020/Scripts/01_SingleOmicsAnalysis"

load(file.path(dataDir, "Transcriptomics/transcriptomics_matrix.Rda"), verbose = TRUE)
load(file.path(dataDir, "Metabolomics/metabolomics_matrix.Rda"), verbose = TRUE)
load(file.path(dataDir, "Cytokines/cytokines_matrix.Rda"), verbose = TRUE)

load(file.path(dataDir, "Transcriptomics/transcriptomics_DE.Rda"), verbose = TRUE)
load(file.path(dataDir, "Metabolomics/metabolomics_DE.Rda"), verbose = TRUE)
load(file.path(dataDir, "Cytokines/cytokines_DE.Rda"), verbose = TRUE)

# Output directories
dir.create(file.path(workDir, "01_clustering"))
dir.create(file.path(workDir, "02_performance"))
dir.create(file.path(workDir, "03_enrichment_analysis"))

# Libraries to use
library(NOISeq)
library(tidyr)
library(ggplot2)
library(plyr)
library(parallel)
library(openxlsx)
library(dplyr)
library(mdgsa)
source(file.path(workDir, "../all_sparse_functions.R"))
source(file.path(workDir, "../auxiliary_functions.R"))


## 1) Clustering #####
# We use DE metabolites + DE cytokines as matrix to clustering analysis in order to reduce the number of PLS models.
metab_sig = metabolomics_DE[metabolomics_DE$adj.P.Val < 0.05,]
metab_sig = t(metabolomics_matrix[rownames(metab_sig),])
cyt_sig = cytokines_DE[cytokines_DE$adj.pval < 0.05,]
cyt_sig = t(cytokines_matrix[rownames(cyt_sig),])

molecules = cbind(cyt_sig, metab_sig)


# Correlation as clustering measure -------------
# We calculate absolute Spearmen correlation. 
# Absolute value because we want high correlation (positive or negative) clusterizing together. 
# The Spearman correlation is less sensitive than the Pearson correlation to strong outliers (human data is very variable with many outlier samples).
cor_matrix = abs(cor(molecules, method = "spearman"))

# Evaluating clustering: best method and best number of clusters by silhouette mean score.  -------------
# pam ----
tested_methods = c("euclidean", "manhattan")
tested_cluster_number = 2:20
cluster_evaluation = do.call(cbind, setNames(lapply(tested_methods, function(i){
  setNames(lapply(tested_cluster_number, function(j){
    out = eval.clusters(matrix = cor_matrix, algorithm = "pam", method = i, clusters = j, corrplot = FALSE)
    out$score
  }), nm = tested_cluster_number)
}), nm = tested_methods))
# plot
data1 = data.frame(tidyr::gather(as.data.frame(cluster_evaluation), key = "method"), tested_cluster_number)
ggplot(data1, aes(fill=method, y=as.numeric(value), x=tested_cluster_number)) + 
  geom_bar(position="dodge", stat="identity") +
  ylab("Average silhouette") + 
  theme(text = element_text(size=20))

# Evaluating clustering: best method and best number of clusters by number of misclassified.  -------------
# pam ----
tested_methods = c("euclidean", "manhattan")
cluster_evaluation = do.call(cbind, setNames(lapply(tested_methods, function(i){
  setNames(lapply(tested_cluster_number, function(j){
    out = eval.clusters(matrix = cor_matrix, algorithm = "pam", method = i, clusters = j, corrplot = FALSE)
    out$bad_classified
  }), nm = tested_cluster_number)
}), nm = tested_methods))

# plot
data2 = data.frame(tidyr::gather(as.data.frame(cluster_evaluation), key = "method"), tested_cluster_number)
ggplot(data2, aes(fill=method, y=as.numeric(value), x=tested_cluster_number)) + 
  geom_bar(position="dodge", stat="identity") +
  ylab("number of misclassified") + 
  theme(text = element_text(size=20))
dev.off()


# Creating clusters: optimal number = 6, best method = euclidean ----
clusters = eval.clusters(matrix = cor_matrix, algorithm = "pam", clusters = 6, method = "euclidean", corrplot = TRUE)$components

# Plot expression profiles
corM = cor(molecules, method = "spearman")

expr.profile.plot(list = clusters, expr.matrix = molecules, corM = corM,
                  path = file.path(workDir, "01_clustering/"))


## 2) PLS performance #####
# Firstly, we define matrices (X, Y) as input:
# As X variable we select differentially expressed genes (statistically different between patients) 
# As Y variable we select cytokines/metabolites of each cluster

# 1.metabolomics
data_met = scale(metab_sig, center = T, scale = F)

# 2.cytokines
data_cyt = scale(cyt_sig, center = T, scale = T)

# 3.transcriptomics
sig_DE_genes = transcriptomics_DE[as.numeric(transcriptomics_DE$adj.P.Val) < 0.05,]
data_rna = transcriptomics_matrix[rownames(sig_DE_genes),]
data_rna = scale(t(data_rna), center = T, scale = F)


# Now, we perform one PLS model per cluster.
X = data_rna # X: transcriptomics (11 patients x 847 genes)
Y = cbind(data_cyt, data_met) # Y: cytokines+metabolomics (11 patients x 35 features)

# This code runs in parallel by every cluster doing a PLS and calculates optimal number of components, R2, Q2 and export all this information as R objects.
mclapply(1:length(clusters), function(i){
  factor = as.factor(metadata$Condition)
  Y = cbind(data_cyt, data_met)
  Y = Y[,clusters[[i]],drop=F] #select variables of each cluster in every loopback
  
  # tune number of components (Ugido's code)
  png(sprintf("02_performance/%02d_%s_r2.png", i, colnames(Y)[1]))
  tune_compR = tune.comp(X, Y, fold = nrow(metadata), factor = factor, ncomp = nrow(metadata), rep = 1, option = "R",
                         title = sprintf("%02d_%s_r2", i, colnames(Y)[1]))
  dev.off()
  png(sprintf("02_performance/%02d_%s_q2.png", i, colnames(Y)[1]))
  tune_compQ = tune.comp(X, Y, fold = nrow(metadata), factor = factor, ncomp = nrow(metadata), rep = 1, option = "Q",
                         title = sprintf("%02d_%s_q2", i, colnames(Y)[1]))
  dev.off()
  
  save(tune_compR, file = sprintf("02_performance/%02d_%s_tunecompR.Rda", i, colnames(Y)[1]))
  save(tune_compQ, file = sprintf("02_performance/%02d_%s_tunecompQ.Rda", i, colnames(Y)[1]))
  
  # Optimal number of components (R2 criteria)
  opt_ncompR = tune_compR$opt.comp
  opt_ncompR_r2 = tune_compR$R2[opt_ncompR]
  opt_ncompR_q2 = tune_compQ$Q2[opt_ncompR]

  # PLS model (Ugidos' code).
  ncomp = opt_ncompR
  pls = pls2(X, Y, a = ncomp, scale = FALSE)

  # Plots
  pdf(sprintf("02_performance/%02d_%s_PLSplots.pdf", i, colnames(Y)[1]))
  par(mfrow=c(2,2))
  score.plot.X(pls, factor)
  loading.plot.X(pls, factor)

  score.plot.Y(pls, factor)
  loading.plot.Y(pls, factor)
  dev.off()

}, mc.cores = 3)

# Making performance table: Number of components, R2, Q2, selected genes
Y = cbind(data_cyt, data_met)
files_var = list.files(path = file.path(workDir, "02_performance/"), pattern = "tunevarQ.Rda")
files_comp = list.files(path = file.path(workDir, "02_performance/"), pattern = "tunecompR.Rda")
files_compQ = list.files(path = file.path(workDir, "02_performance/"), pattern = "tunecompQ.Rda") # fix to extract Q2

performance_table = do.call(rbind, lapply(1:length(files_comp), function(i){

    # sPLS model
    load(file = file.path(workDir, "02_performance/", files_comp[i]))
    load(file = file.path(workDir, "02_performance/", files_var[i]))
    load(file = file.path(workDir, "02_performance/", files_compQ[i]))

    opt_ncomp = tune_compR$opt.comp
    opt_nvar = tune.varQ$var.opt.x
    
    spls = spls_fx(X, Y[,clusters[[i]],drop=F], n=opt_ncomp, kx = opt_nvar)
    selected_genes = rownames(spls$W[rowSums(abs(spls$W)) != 0,])
    non_selected_genes = rownames(spls$W[rowSums(abs(spls$W)) == 0,])

    opt_var_value = length(selected_genes)

    # We are extracting the list selected and non-selected genes by PLS in each cluster for further steps
    # write.table(sort(selected_genes), file = file.path(workDir, "02_performance", sub("_tunecompQ.Rda", replacement = "_genesSELECTED.txt", x = files_comp[i])),
    #             quote = F, row.names = F, col.names = F)
    # write.table(sort(non_selected_genes), file = file.path(workDir, "02_performance", sub("_tunecompQ.Rda", replacement = "_genesNONselected.txt", x = files_comp[i])),
    #             quote = F, row.names = F, col.names = F)
    
    tabla2 = data.frame("Number of components" = opt_ncomp,
                        "R2" = round(tune_compR$R2[opt_ncomp], 3),
                        "Q2" = round(tune_compQ$Q2[opt_ncomp], 3),
                        "Selected genes by Q2" = opt_var_value,
                        #"%lnc" = num_lnc,
                        #"%TF" = num_TF,
                        row.names = files_comp[i], check.names = FALSE)
    return(tabla2)
    }))

save(performance_table, file = file.path(workDir, "02_performance/performance_table.Rda"))


## 3) Biological interpretation #####
# Enrichment analysis of every group of genes selected in each PLS model: 
# Assuming every cluster of metabolites/interleukins represents one biological function, we want to infer that function by PLS selected genes for that cluster.

# annotation: GO ----
load("/home/teresa/Desktop/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_GOannotation.Rda", verbose=T)
myannot = myannot[myannot$hgnc_symbol %in% colnames(data_rna),]
anot2 = as.data.frame(myannot %>% group_by(go_id) %>% filter(n() >= 10))
anot2 = as.data.frame(anot2 %>% group_by(go_id) %>% filter(n() <= 500))

anot = anot2[anot2$namespace_1003 %in% c("biological_process"),]
GO_bp = anot[,c("hgnc_symbol", "name_1006", "go_id")]

annotation = list(GO_bp = GO_bp)

## 3.1) mdgsa with PC1 loadings ----
Y = cbind(data_cyt, data_met)
files = list.files(path = file.path(workDir, "02_performance"), pattern = "tunecompQ.Rda")
annot = annotMat2list (GO_bp[,c(1,2)]) #Already filtered

png(file.path(workDir, "../../MHEarticle/Figures/S3_Figure.png"), width = 674*3, height = 552*3)
par(mfrow=c(3,3))
EA = lapply(1:length(clusters), function(i){#c(1,2,3,4,6), function(i){
    
    # Extracting significative genes
    pls = pls2(X, Y[,clusters[[i]],drop=F], a = performance_table[i,1], scale = FALSE)
    
    sig_genes = pls$W[,1,drop=F] # 1st component X loadings
    components = c(1,2)
    
    # sig_genes = sig_genes[!rownames(sig_genes) %in% colnames(lnc),,drop=F] # Keep out lncRNAs
    par(mar=c(6,6,6,6)+.1) 
    withoutMHE = sum(metadata$Condition == "withoutMHE")
    if (mean(pls$T[1:withoutMHE, 1]) < 0){ #If the mean of group nMHE is negative, we keep loadings with de same sign
      pls$T = pls$T
    } else { #If the mean of group nMHE es positive, we need to change the sign
      pls$T = - pls$T
    }
    
    plot(pls$T[,components[1]], pls$T[,components[2]], pch = 19, cex = 8, cex.axis = 3, cex.lab = 4, cex.main = 4,
         col = c("gold","mediumblue")[factor(metadata$Condition)], 
         main = paste("Module", i),
         xlab = paste("Comp", components[1]), ylab = paste("Comp", components[2]),
         xlim = c(-20,20), ylim = c(-20,20))
    # text(x = pls$T[,components[1]], y = pls$T[,components[2]]-1, labels = names(pls$T[,1]), pos = 1, cex = 4)
    if(nrow(sig_genes) > 1) {
      # We use loadings as p.values. It is not necessary to use pval2index() function because loadings it is already the rindex we need.
      # The logistic model works well even if rindex is not a value between 0-1 (because when log=T, the pval achieves values bigger than 1).
      # rindex = pval2index (pval = 1 + abs(sig_genes[,1,drop=F]), sign = sign(sig_genes[,1,drop=F]), log = FALSE)

      # Comparison MHE vs nMHE: because sometimes PLS outputs loadings in the oposite direction, we need to check the sign
      withoutMHE = sum(metadata$Condition == "withoutMHE")
      if (mean(pls$T[1:withoutMHE, 1]) < 0){ #If the mean of group nMHE is negative, we keep loadings with de same sign
        rindex = sig_genes
      } else { #If the mean of group nMHE es positive, we need to change the sign
        rindex = - sig_genes
      }

      # Transformation of the ranking index so that its distribution is suitable as independent variable of the logistic regression model.
      rindex = indexTransform (rindex)

      #Now everything is ready to carry out the univariate gene set analysis.
      res.uv = uvGsa_withGeneNames(rindex, annot, fulltable = T, family = gaussian())#, control = list(maxit=30))
      res.uv = res.uv[order(as.numeric(res.uv$pval)),]

    } else {
      res.uv = data.frame("N"=NA, "lor"=NA, "pval"=NA, "padj"=NA, "sd"=NA, "t"=NA, "conv"=NA)
    }

    write.table(res.uv, file=file.path(workDir, "03_enrichment_analysis/mdgsa",
                                       paste0("Cluster", sprintf("%02d", i), "_", clusters[[i]][1], "_GO")),
                                       sep = "\t", col.names=NA, quote = F)

  })
dev.off()

## Summary and representation ----
directories = "mdgsa"

# Datasets per sheet
list_of_datasets = setNames(lapply(directories, function(i){
  
  files = list.files(path = file.path(workDir, "03_enrichment_analysis/", i), 
                      include.dirs = F, pattern = "_GO$")

  dfs = do.call(rbind, lapply(files, function(x){
    df = read.delim(file.path(paste0("03_enrichment_analysis/", i), x), header = F, stringsAsFactors = F)
    header = df[1,]
    tmp_df = df[-1,]
    df = tmp_df[which(as.numeric(tmp_df$V3) < 0.05),] #V3: pvalue, V4: adjPval
    df = rbind(header, df)
    df = rbind(c(x,rep(NA,ncol(df)-1)), df)
    df[c(nrow(df)+1, nrow(df)+2),] = NA
    df
  }))
  
  return(dfs)
  
}), nm = directories)

# Write excel workbook
write.xlsx(list_of_datasets, file = file.path(workDir, "03_enrichment_analysis/mdgsa/EA_mdgsa_GO_summary__X1090DE_pam_6clusters.xlsx"))


# Plots ----
myplots = setNames(lapply(directories, function(i){
  
  files = list.files(path = paste0("03_enrichment_analysis/", i), include.dirs = F, pattern = "GO$")
  
  dfs = do.call(rbind, lapply(files, function(x){
    df = read.delim(file.path("03_enrichment_analysis/", i, x), header = T)
    df = df[which(as.numeric(df$pval) < 0.05),]
    if (length(df$N) < 1) {
      a = data.frame("df.N"=0, x)
    } else {
      a = data.frame(df$N, x)
    }
    a
  }))
  ## barplot----
  # png(filename = file.path(path, i, "barplot_mdgsa.png"), width = 7180, height = 3340, res = 500)
  barplot(dfs$df.N, 
          axes = F,
          ylab = "Number of annotated genes",
          col = pals::glasbey()[-c(1,2)][factor(dfs$x)])
  axis(2, at = seq(0, max(dfs$df.N)+1, by = 1))
  # abline(h = 5, col = "red")
  # dev.off()
  
  # Mascleta plot ----
  # First input to cytoscape: 1 hub per enrichment linked to each significative pathway
  dfs_partA = do.call(rbind, lapply(files, function(x){
    df = read.delim(file.path(workDir, "03_enrichment_analysis/", i, x), header = F, stringsAsFactors = F)
    # Keeping columns: term, annotTest, pval and annotTest_names
    df = df[,c(1,8,3,9)]
    # Adding a column to link each pathway to the same hub
    cluster = substr(gsub(paste0(i,"_"), "",x), start = 1, stop = 9)
    df$V5 = cluster
    # Changing header with hub parameters: "term = hub name", "annotTest = node size = 0", "pval=1", "annotTest_names=NULL", "link to itself = NULL"
    df = df[-1,]
    df = df[as.numeric(df[,3]) < 0.05,]
    df = rbind(c(cluster, 0, 1, " ", cluster), df)
    
    colnames(df) = paste0("V", seq(1:5))
    df
  }))
  # Second input to cytoscape: shared pathways between hubs (pie charts)
  # create the dummies
  dfs_partB = tidyr::spread(dfs_partA[,c(1,2,5)], key = V5, value = V2)
  dfs_partB$sum = rowSums(do.call(cbind, lapply(dfs_partB[,-1], as.numeric)), na.rm = T)
  # end
  return(list(dfs_partA, dfs_partB))

}), nm = directories)

# Write excel workbook
write.xlsx(myplots$mdgsa[[1]], file = file.path("03_enrichment_analysis/", directories, "EA_mdgsa_GO_mascletaplot_partA_newbp.xlsx"), colnames = F)
write.xlsx(myplots$mdgsa[[2]], file = file.path("03_enrichment_analysis/", directories, "EA_mdgsa_GO_mascletaplot_partB_newbp.xlsx"), colnames = T)

