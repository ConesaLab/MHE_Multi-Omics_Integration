#######################################################################
#######             METABOLOMICS ANALYSIS MHE DATA              #######
#######################################################################

## By Teresa Rubio
# Raw data are available at supplementary material: S7_Table.xlsx
# Metabolomics platform: 'AbsoluteIDQ® p180 Kit' (Biocrates)
# List of Metabolites: https://biocrates.com/wp-content/uploads/2020/02/Biocrates_p180_metabolites.pdf
# Metabolomics service: 'Protein and Metabolite Biomarkers Unit, Istituto di Ricerche Farmacologiche Mario Negri IRCCS (Milano, Italy)'

  
# Work directory
workDir = "~/Desktop/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Metabolomics"
setwd(workDir)
options(stringsAsFactors = FALSE)

# Libraries to use
library(readxl)
library(vsn)
library(limma)
library(NOISeq)
source(file.path(workDir, "../../auxiliary_functions.R"))


## Step 1:  Raw data quality control and normalization ###### 

  # 1. Data matrix and Metadata ----
  S6_Table <- read_excel("~/Desktop/MHERubio2020/MHEarticle/Final_figures/S6_Table.xlsx")
    
  metadata = S6_Table[, c("SampleName", "Condition")]
  raw_data = t(S6_Table[, -c(1,2)])
  colnames(raw_data) = metadata$SampleName
  
  # 2. Quality control
  # PCA (NOISeq) ----
  myfactors = data.frame(Sample = metadata$SampleName, Type = metadata$Condition)
  exprdata = readData(data=raw_data, factors=myfactors)                   
  myPCA = dat(exprdata, type = "PCA", logtransf = F, norm = T) #Apply log-transformation to improve plots
    
  plot.scores.gg(myPCA, comp = c(1,2), title = "Score plot - Raw data")
  plot.loadings.gg(myPCA, comp = c(1,2), title = "Loading plot - Raw data")
  
  # Boxplot ----
  mycolors = c("salmon", "lightblue")[as.factor(metadata$Condition)]
  boxplot(raw_data+1, las=2, log="y", col=mycolors, names = metadata$SampleName,
          main = "Metabolomics before normalization", ylab = "log scale")
  legend("topright", y =c("withoutMHE", "whithMHE"), pch = 15, col = c("lightblue", "salmon"))
    
  # 3. Normalization  ----
  vsn = normalizeVSN(raw_data) # --> log2-normalized data
  # output
  metabolomics_matrix = vsn
  save(metabolomics_matrix, file = file.path(workDir, "metabolomics_matrix.Rda"))
  
  # 4. Quality control after normalization ----
  # PCA (NOISeq) ----
  myfactors = data.frame(Sample = metadata$SampleName, Type = metadata$Condition)
  exprdata = readData(data=vsn, factors=myfactors)                   
  myPCA = dat(exprdata, type = "PCA", logtransf = T, norm = T) #Already log-transformed and normalized data
  
  plot.scores.gg(myPCA, comp = c(1,2), title = "Score plot - VSN normalized data")
  plot.loadings.gg(myPCA, comp = c(1,2), title = "Loading plot - VSN normalized data")
  
  # Boxplot ----
  mycolors = c("salmon", "lightblue")[as.factor(metadata$Condition)]
  boxplot(vsn, las=2, col=mycolors, names = metadata$SampleName,
          main = "Metabolomics after normalization")
  legend("topright", y =c("withoutMHE", "whithMHE"), pch = 15, col = c("lightblue", "salmon"))
  
  
## Step 2: Metabolomis analysis ###### 
    
  # 1. Differential expression ----
  # Build the design matrix for the linear modelling function.
  f = factor(metadata$Condition, levels = unique(metadata$Condition))
  design = model.matrix(~0 + f)
  colnames(design) = c("without_MHE", "with_MHE")
  
  # Apply lmFit to the normalized values.
  fit = lmFit(vsn, design)
  
  # Create a contrast matrix.
  contrast.matrix = makeContrasts("with_MHE-without_MHE", levels = design)
  
  # Apply this contrast matrix to the modeled data and compute statistics for the data.
  fit = contrasts.fit(fit, contrast.matrix)
  fit = eBayes(fit)
  
  metabolomics_DE = topTable(fit, coef = 1, number = Inf, adjust.method = "none")
  
  save(metabolomics_DE, file = file.path(workDir, "metabolomics_DE.Rda"))
    
          
  # 2. Enrichment analysis: PaintOmics ----
  ## Data file ----
    # scale metabolite values due magnitude differences
    scaled_metabolomics = scale(raw_data, center = T, scale = T)
    # mean expression values for each group of patients
    mean_table = aggregate(t(scaled_metabolomics), by = list(metadata$Condition), FUN = mean)
    rownames(mean_table) = mean_table$Group.1
    mean_table = mean_table[,-1]
    mean_table = t(mean_table)
    # calculate mean difference between groups of patients
    diff = matrix(ncol = 2, nrow = nrow(mean_table))
    colnames(diff) = c("feature_names", "withMHE.vs.withoutNMHE")
    diff[,"withMHE.vs.withoutNMHE"] = mean_table[,"withMHE"] - mean_table[,"withoutMHE"]
    diff[,"feature_names"] = rownames(mean_table)
    # export file with PaintOmics format
    write.table(diff, file=file.path(workDir, "metabolomics_Paintomics_datafile.csv"), sep="\t", row.names = F)
    
  ## Relevant features file ----
    # we select differentially expressed genes from limma test
    significant_metabolites <- metabolomics_DE[which(metabolomics_DE$adj.P.Val < 0.05),]
    relevant_features <- data.frame("feature_names" = rownames(significant_metabolites))
    write.table(relevant_features, file=file.path(workDir, "metabolomics_Paintomics_relevantfeatures.csv"), sep="\t", row.names = F)
    
  
