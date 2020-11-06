These are the instruction to reproduce the results of the paper called 
**“Multi-omic analysis unveils biological pathways in peripheral immune system associated to minimal hepatic encephalopathy appearance in cirrhotic patients”**

*Teresa Rubio, Vicente Felipo, Sonia Tarazona, Roberta Pastorelli, Desamparados Escudero-García, Joan Tosca, Amparo Urios, Ana Conesa, Carmina Montoliu.*


**Scripts**

	- **01_SingleOmicsAnalsis**: 
	  - transcriptomics (R script for differential expression + enrichment analysis): output + R script
	  - metabolomics (R script for differential expression analysis): output + R script
	  - cytokines (R script for differential expression analysis): output + R script

	- **02_Multi-Omics_Integration**:
	  - 01_clustering (output of metabolites and cytokines clustering analysis)
	  - 02_performance (output of PLS analysis and model performance)
	  - 03_enrichment_analysis (output of mdgsa analysis with genes ranked by PLS loadings)
	  - MultiOmicsAnalysis_script.R (R script to produce the 3 above outputs)
	  
	- **revision**: scripts to answer reviewers questions
	
	- **all_sparse_functions.R** (R script with homemade PLS additional functions)
	
	- **auxiliary_functions.R** (R script with homemade auxiliary functions)
	
