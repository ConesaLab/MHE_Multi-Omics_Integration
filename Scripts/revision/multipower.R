#######################################################################
#######               Multipower ANALYSIS MHE DATA              #######
#######################################################################

## By Teresa Rubio

# Work directory
workDir = "~/ClusterCIPF/MHERubio2020/Scripts/revision"
setwd(workDir)
options(stringsAsFactors = FALSE)

# Libraries to use
library(MultiPower)

# 1. Data matrix and Metadata ----
load("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_matrix.Rda")
load("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Metabolomics/metabolomics_matrix.Rda")
load("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Cytokines/cytokines_matrix.Rda")

metadata = read.delim("~/ClusterCIPF/MHERubio2020/Scripts/01_SingleOmicsAnalysis/Transcriptomics/transcriptomics_patients.txt")
metadata = metadata[,c("SampleName", "Condition")]


# Calculating statistical power by hand ---
# function power.t.test
data = cytokines_matrix


# 1. Calculate delta as difference between means
geneMean = t(apply(data, 1, tapply, metadata$Condition, mean))
delta = abs(geneMean[,"withMHE"] - geneMean[,"withoutMHE"])
delta = sort(delta, decreasing = T)

delta1 = median(delta[1:(length(delta) * 0.01)])
delta5 = median(delta[(length(delta) * 0.01):(length(delta) * 0.05)])
delta10 =  median(delta[(length(delta) * 0.05):(length(delta) * 0.10)])

delta25 = median(delta[(length(delta) * 0.10):(length(delta) * 0.25)])
delta50 = median(delta[(length(delta) * 0.25):(length(delta) * 0.50)])
delta75 = median(delta[(length(delta) * 0.50):(length(delta) * 0.75)])
delta100 = median(delta[(length(delta) * 0.75):(length(delta))])


# 2. Calculate s as pooled standard deviation
# Standard deviation per group
groups = metadata$Condition
sdPerGroup = t(apply(data, 1, tapply, INDEX = groups, sd, na.rm = TRUE))

# Pooled Standard Deviation
nGroup = table(groups)
SDpooled = sqrt((nGroup[1]*sdPerGroup[,1]^2 + nGroup[2]*sdPerGroup[,2]^2)/(sum(nGroup)-2))
SDpooled = SDpooled[names(delta)]


# 3. Percentiles
quantiles = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)

deltas = stats::quantile(delta, probs = quantiles)
eses = stats::quantile(SDpooled, probs = quantiles)

data.frame ("delta" = deltas, "eses" = eses)

names(deltas) = paste("delta", names(deltas))
names(eses) = paste("s", names(eses))


# Calculate Cohen's d
d = do.call(cbind, lapply(deltas, function(d){
  do.call(rbind, lapply(eses, function(s){
    d / s
  }))
}))

# Calculate power
power = do.call(cbind, lapply(deltas, function(d){
  do.call(rbind, lapply(eses, function(s){
    potencia = power.t.test(n = 6, delta = d, sd = s, sig.level = 0.05,
                            type = "two.sample", alternative = "two.sided")
    potencia$power
  }))
}))

# colnames
colnames(d) = names(deltas)
colnames(power) = names(deltas)

d
power

