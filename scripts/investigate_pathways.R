library(data.table)

data = fread("results/pathway_analysis/c7_immunology/intNMF/component2_pathways.tsv")
data

colnames(data)
dim(data)
data[, 1:7]
