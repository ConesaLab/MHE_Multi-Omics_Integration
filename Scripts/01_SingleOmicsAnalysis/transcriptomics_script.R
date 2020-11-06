#######################################################################
#######           TRANSCRIPTOMICS ANALYSIS MHE DATA             #######
#######################################################################

## By Teresa Rubio
# Raw data are available at GEO: GSE149741
# Annotation available at GEO: GPL20844 (https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL20nnn/GPL20844/suppl/GPL20844_072363_D_GEO_20150612.txt.gz)
# Hybridization platform: 'SurePrint G3 Human Gene Expression Microarrays v3 8X60K [feature version]' (Agilent p/n G4858A-072363)
# Sequentiation service: 'Servicio de Genomica y Genetica Translacional, CIPF (Valencia, Spain)'

  
# Working directory
workDir = "~/Desktop/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics"
setwd(workDir)
options(stringsAsFactors = FALSE)

# Libraries to use
library(GEOquery)
library(vsn)
library(readr)
library(limma)
library(NOISeq)
library(biomaRt)
library(mdgsa)
source(file.path(workDir, "../../auxiliary_functions.R"))


## Step 1:  Raw data to expression matrix ###### 

  # 1. Metadata ----
    # The tab-delimited file that contains both the names of your Agilent Feature Extraction 
    # raw data text files and the corresponding sample information.
    targets = readTargets(file.path(workDir, "transcriptomics_patients.txt"), sep="\t", row.names = "SampleName")


  # 2. Download rawdata from GSE149741  ----
    # filePaths = getGEOSuppFiles("GSE149741")
    # # Define data directory
    # dataDir = gsub(pattern = "[^/]+$*", replacement = "", x = rownames(filePaths)[2])
    # # Extract compressed files
    # untar(tarfile = rownames(filePaths)[2], exdir = dataDir)
    # lapply(file.path(dataDir, list.files(dataDir, pattern = ".gz")), FUN = gunzip, remove = T)
    # # Use LIMMA's read.maimages function to load raw data into an RGList object. Remember to set the path to the 
    # # location where your Agilent FE text files are stored.
    # raw_data = read.maimages(targets, path=file.path(dataDir), source = "agilent", green.only=TRUE)
    raw_data = read.maimages(targets, path=file.path("GSE149741/"), source = "agilent", green.only=TRUE)
    colnames(raw_data$E) = targets$SampleName
    rownames(raw_data$E) = raw_data$genes$ProbeName

    
  # 3. Normalization  ----
    # vsn algorithm performs background correction + vsn normalization simultaneously --> log2-normalized
    vsn = normalizeVSN(raw_data)
    
    
  # 4. Filtering low expressed genes  ----
    # To get an idea of how bright expression probes should be, 
    # we compute the 80% percentile of the negative control probes on each array. 
    # We keep probes that are at least 10% brighter than the negative controls 
    # on at least 5 arrays (because there are 5 patients per group)
    neg95 = apply(vsn$E[vsn$genes$ControlType==-1,], 2, function(x) quantile(x, p = 0.8))
    cutoff = matrix(1.1 * neg95, nrow(vsn), ncol(vsn), byrow = TRUE)
    isexpr = rowSums(vsn$E > cutoff) >= 5

    vsn_filtered <- vsn[vsn$genes$ControlType==0 & isexpr,]
    
    
  # 5. Averaging replicate spots  ----
    agilent.db = read_delim(file.path(workDir, "GPL20844_072363_D_GEO_20150612.txt"), 
                            "\t", escape_double = FALSE, comment = "!", trim_ws = TRUE, skip = 42)
    # remove control spots. 
    # In Agilent annotation, "CONTROL_TYPE = pos" are positive control spots, "CONTROL_TYPE = neg" are negative control spots and "CONTROL_TYPE = FALSE" are probes.
    agilent.db = agilent.db[agilent.db$CONTROL_TYPE==FALSE,]

    # this GEO annotation has errors due to Excel formatting that were transported into the platform record
    unique(agilent.db[order(agilent.db$GENE_SYMBOL),]$GENE_SYMBOL[1:55])
    
    # we fix those gene symbols manually
    agilent.db$GENE_SYMBOL = gsub("^1-Dec", "DELEC1", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^1-Mar", "MARCHF1", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^1-Sep", "SEPTIN1", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^10-Mar", "MARCHF10", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^10-Sep", "SEPTIN10", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^11-Mar", "MARCHF11", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^11-Sep", "SEPTIN11", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^12-Sep", "SEPTIN12", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^14-Sep", "SEPTIN14", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^15-Sep", "SEPTIN15", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^2-Mar", "MARCHF2", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^2-Sep", "SEPTIN2", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^3-Mar", "MARCHF3", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^3-Sep", "SEPTIN3", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^4-Mar", "MARCHF4", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^4-Sep", "SEPTIN4", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^5-Mar", "MARCHF5", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^5-Sep", "SEPTIN5", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^6-Mar", "MARCHF6", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^6-Sep", "SEPTIN6", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^7-Mar", "MARCHF7", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^7-Sep", "SEPTIN7", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^8-Mar", "MARCHF8", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^8-Sep", "SEPTIN8", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^9-Mar", "MARCHF9", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    agilent.db$GENE_SYMBOL = gsub("^9-Sep", "SEPTIN9", agilent.db$GENE_SYMBOL, ignore.case = FALSE)
    
    # merge expression matrix with annotation
    transcriptomics_matrix = vsn_filtered$E
    transcriptomics_matrix = cbind(rownames(transcriptomics_matrix), transcriptomics_matrix)
    colnames(transcriptomics_matrix)[1] = "NAME"
    transcriptomics_matrix = merge(transcriptomics_matrix, agilent.db[,c("NAME", "GENE_SYMBOL")], by = "NAME")
    # each gene may be represented by one or more probes, replicate probes are replaced with their average
    transcriptomics_matrix = avereps(transcriptomics_matrix, ID = transcriptomics_matrix$GENE_SYMBOL)
    transcriptomics_matrix = transcriptomics_matrix[-c(which(is.na(transcriptomics_matrix[,c("GENE_SYMBOL")]))),] #remove empty names
    
    rownames(transcriptomics_matrix) = transcriptomics_matrix[,c("GENE_SYMBOL")]
    transcriptomics_matrix = transcriptomics_matrix[,-c(1,ncol(transcriptomics_matrix))]
    transcriptomics_matrix = apply(transcriptomics_matrix, c(1,2), as.numeric)

    save(transcriptomics_matrix, file = file.path(workDir, "transcriptomics_matrix.Rda"))
    
    
## Step 2: Transcriptomics analysis ###### 
    
  # 1. PCA (NOISeq) ----
    myfactors = data.frame(Sample = targets$SampleName, Type = targets$Condition)
    exprdata = readData(data=transcriptomics_matrix, factors=myfactors)                   
    myPCA = dat(exprdata, type = "PCA", logtransf = T, norm = T) #Already log-transformed and normalized data
  
    plot.scores.gg(myPCA, comp = c(1,2), title = "Score plot - VSN normalized data")
    plot.loadings.gg(myPCA, comp = c(1,2), title = "Loading plot - VSN normalized data")
  
    
  # 2. Differential expression ----
    f = factor(targets$Condition, levels = unique(targets$Condition))
    design = model.matrix(~0 + f)
    colnames(design) = c("without_MHE", "with_MHE")
    
    # Calculate weights
    arrayw = arrayWeights(transcriptomics_matrix)
    
    # Apply the intensity values to lmFit.
    fit = lmFit(transcriptomics_matrix, design, weights = arrayw)
        
    # Create a contrast matrix.
    contrast.matrix = makeContrasts("with_MHE-without_MHE", levels = design)
        
    # Apply this contrast matrix to the modeled data and compute statistics for the data.
    fit = contrasts.fit(fit, contrast.matrix)
    fit = eBayes(fit, trend = F, robust = F)
    
    transcriptomics_DE = topTable(fit, coef = 1, number = Inf, adjust.method = "none")
    
    save(transcriptomics_DE, file = file.path(workDir, "transcriptomics_DE.Rda"))
    
          
    # 3. Enrichment analysis ----
    # database used: GO (biological process)
    biomartHuman = useMart("ensembl", dataset = "hsapiens_gene_ensembl", version="Ensembl Genes 100") #host = "apr2020.archive.ensembl.org"
    myannot = getBM(attributes = c("hgnc_symbol","go_id", "name_1006", "namespace_1003", "definition_1006"),
                    filters = "hgnc_symbol", 
                    values = rownames(transcriptomics_matrix),
                    mart = biomartHuman)
    save(myannot, file = file.path(workDir, "transcriptomics_GOannotation.Rda"))
    
    myannot = myannot[myannot$namespace_1003 == "biological_process",]
      
    # enrichment analysis method: mdgsa
    rindex = pval2index (pval = fit$p.value, sign = fit$t)
    rindex = indexTransform (rindex)
           
    annot = annotMat2list (myannot[,c(1,3)])
    annot = annotFilter (annot, rindex, minBlockSize = 10, maxBlockSize = 500, verbose = TRUE)
           
    transcriptomics_EA = uvGsa_withGeneNames(rindex, annot)
    
    save(transcriptomics_EA, file = file.path(workDir, "transcriptomics_EA.Rda"))
    
    # input to REVIGO web (to summarize EA output)
    transcriptomics_EA$name_1006 = rownames(transcriptomics_EA)
    transcriptomics_EA = merge(unique(myannot[,2:3]), transcriptomics_EA)
    transcriptomics_EA = transcriptomics_EA[order(transcriptomics_EA$pval),]
         
    sig_up = transcriptomics_EA[which(transcriptomics_EA$lor > 0 & transcriptomics_EA$pval < 0.05),]
    sig_down = transcriptomics_EA[which(transcriptomics_EA$lor < 0 & transcriptomics_EA$pval < 0.05),]
    
    options(scipen=10) # quit scientific notation
    write.csv(sig_up[, c("go_id","pval")], file = file.path(workDir, "transcriptomics_revigoUP_input.csv"), row.names = F)
    write.csv(sig_down[, c("go_id","pval")], file = file.path(workDir, "transcriptomics_revigoDOWN_input.csv"), row.names = F)
    options(scipen=0)  # restore the default
    
    # REVIGO website parameters:
    # cutoff <- "0.5" #Allowed values: "0.90" (large) "0.70" (normal) "0.50" (small) "0.40" (tiny) 
    # organism <- "whole UniProt" #Allowed values: See organism.list below
    # isPValue <- "yes" #Allowed values: "yes"  "no"
    # whatIsBetter <- "higher" #Allowed values: "higher" "lower" "absolute" "abs_log"
    # measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"
    # domain <- "BP" #Allowed values: "BP" "CC" "MF"
    
    
# EXTRA: PaintOmics input files to help with interpretation ###### 
  ## Data file ----
    # mean expression values for each group of patients
    mean_table = aggregate(t(transcriptomics_matrix), by = list(targets$Condition), FUN = mean)
    rownames(mean_table) = mean_table$Group.1
    mean_table = mean_table[,-1]
    mean_table = t(mean_table)
    # calculate mean difference between groups of patients
    diff = matrix(ncol = 2, nrow = nrow(mean_table))
    colnames(diff) = c("feature_names", "withMHE.vs.withoutNMHE")
    diff[,"withMHE.vs.withoutNMHE"] = mean_table[,"withMHE"] - mean_table[,"withoutMHE"]
    diff[,"feature_names"] = rownames(mean_table)
    # export file with PaintOmics format
    write.table(diff, file=file.path(workDir, "transcriptomics_Paintomics_datafile.csv"), sep="\t", row.names = F)
    
  ## Relevant features file ----
    # we select differentially expressed genes from limma test
    significant_genes <- transcriptomics_DE[which(transcriptomics_DE$adj.P.Val < 0.05),]
    relevant_features <- data.frame("feature_names" = rownames(significant_genes))
    write.table(relevant_features, file=file.path(workDir, "transcriptomics_Paintomics_relevantfeatures.csv"), sep="\t", row.names = F)
