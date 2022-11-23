library("ggplot2")
library("reshape2")
library("stringr")
library("ggpubr")
library("gridExtra")
library("dplyr")
library("ComplexHeatmap")
library("circlize")
library('biomaRt')

args <- commandArgs()
args = c("data/cmi_pb/cytof.pma.day0.proc.tsv", 
         "data/cmi_pb/olink.pma.day0.proc.tsv",
         "data/cmi_pb/rnaseq.pma.day0.proc.tsv",
         "results/colab_anna_max/cytof.loadings.tsv",
         "results/colab_anna_max/olink.loadings.tsv",
         "results/colab_anna_max/rnaseq.loadings.tsv",
         "results/colab_anna_max/factors_mcia.clean.tsv",
         "results/colab_anna_max/factors_v_features/")
cytof_cmidb_fn = args[1]
olink_cmidb_fn = args[2]
rnaseq_cmidb_fn = args[3]
cytof_loading_fn = args[4]
olink_loading_fn = args[5]
rnaseq_loading_fn = args[6]
factors_fn = args[7]
results_folder = args[8]

# data to rename cytof columns
cytof_dict = read.table("results/cytof_dict.tsv", row.names = 1)

# creating the results folder
dir.create(results_folder, showWarnings = F)

# loading the factor data
factors = read.table(factors_fn, row.names = 1, header=1, sep='\t')

# loading the features, loadings and metagenes
features = list()
loadings = list()
metagenes = list()

for (assay in c('cytof', 'olink', 'rnaseq')){
    
    # load the original feature information
    orig_fn = paste0("data/cmi_pb/", assay, ".pma.day0.proc.tsv")
    orig = read.table(orig_fn, header=T, sep='\t', row.names = 1)
    orig = orig[match(rownames(factors), rownames(orig)), ]
    features[[assay]] = orig
    
    # load the loading values
    loading_fn = paste0("results/colab_anna_max/", assay, ".loadings.tsv")
    loadings = read.table(loading_fn, header=T, quote = "", sep='\t', row.names = 1)
    
    if (assay == 'cytof'){
        reorder = match(rownames(loadings), rownames(cytof_dict))
        rownames(loadings) = cytof_dict[reorder,]
    }
    else if (assay == 'rnaseq'){
        rownames(loadings) = str_replace(rownames(loadings), '\\.[0-9]+$', '')
    }
    
    loadings[[assay]] = loadings
    
    # calculate the metagene values
    shared_features = intersect(colnames(orig), rownames(loadings))
    shared_orig = orig[, shared_features]
    shared_orig = orig[, order(colnames(shared_orig))]
    shared_loadings = loadings[shared_features, ]
    shared_loadings = loadings[order(rownames(shared_loadings)), ]
    calc_metagene = as.matrix(shared_orig) %*% as.matrix(shared_loadings)
    metagenes[[assay]] = calc_metagene
    
}

# loading the biomart
ensembl = useMart("ensembl", host="useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", ensembl)

##################################
##### Preprocessing the data #####
##################################
# change the loading ensembl id to hgnc name
genes = colnames(features[['rnaseq']])
hgnc_list <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   mart=mart)

hgnc_list = as.data.frame(hgnc_list)
hgnc_list = hgnc_list[(hgnc_list[2] != ""), ]

print("obtaining hgnc name")
print("renaming now")
reorder = match(colnames(features[['rnaseq']]), 
                         hgnc_list[['ensembl_gene_id']])
colnames(features[['rnaseq']]) = hgnc_list[reorder, 'hgnc_symbol']
features[['rnaseq']] = features[['rnaseq']][, !(is.na(colnames(features[['rnaseq']])))]

###########################################################
#### Plot the correlation between factors and features ####
###########################################################

# Function: correlates factors versus feature data
make_heatmap <- function(feat,
                         feat_data,
                         factors,
                         clust_cols=T,
                         clust_rows=T, 
                         vis_fill_vals=F){
  
  # set the cell func
  if (vis_fill_vals == T){
    cell_func <- function(j, i, x, y, width, height, fill){
      text = sprintf("%.2f", factor_corrs[i, j])
      grid.text(text, x, y, gp = gpar(fontsize = 5))
    }
  }
  else {
    cell_func <- NULL
  }

  # correlate the factors and task
  print("# correlate the factors and features")
  factor_corrs = cor(factors, feat_data)
  factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
  color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # make a heatmap of the correlations
  print("# make a heatmap of the correlations")
  title = sprintf('factors versus %s', toupper(feat))
  show_cols = T
  if (feat %in% c('rnaseq', 'olink')){
    show_cols = F
  }
  
  p = Heatmap(factor_corrs, 
          name = "Pearson's R", 
          column_title = title,
          row_title = "Factors",
          row_names_gp = gpar(fontsize = 7),
          col = color_func,
          cluster_rows = clust_rows,
          show_column_names = show_cols,
          cluster_columns = clust_cols,
          cell_fun = cell_func,
          heatmap_width = unit(5.5, "in"),
          heatmap_height = unit(6, "in"))
  
  return(p)
}

for (feat_name in names(features)){
  
  # save the heatmaps to separate pdf pages
  fn = paste0(results_folder, feat_name, '.pdf')
  pdf(fn, onefile = TRUE)
  
  # correlate the factors and task
  print("# correlate the factors and features")
  factor_corrs = cor(factors, features[[feat_name]])
  print(dim(factor_corrs))
  factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
  print(dim(factor_corrs))
  
  color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # make a heatmap of the correlations
  print("# make a heatmap of the correlations")
  title = sprintf('%s factors versus %s', 'mcia', toupper(feat_name))
  
  if (feat_name %in% c("rnaseq", "olink")){
    vis_fill_vals = F
  }
  else{
      vis_fill_vals = T
  }
  p = make_heatmap(feat_name,
           features[[feat_name]],
           factors,
           method=feat_name, 
           clust_row = F,
           clust_cols = clust_cols,
           vis_fill_vals = vis_fill_vals)
  print(p)
  
  # close pdf
  dev.off()
}

#####################################################
####### Plot RNA-seq genes with most variance########
#####################################################

## filtering based on the variances of gene expression
## filtering is done using a percentile min and max
#rnaseq_variance = rnaseq %>% summarise_all(var)
#rnaseq_variance =  t(rnaseq_variance)

## Returns a vector for easy filtering
## filter based on the percentile using pmin and pmax
#filter_by_percentile <- function(vec, pmin, pmax){
#  var_func = ecdf(rnaseq_variance[, 1])
#  bools = list()
#  for (i in seq(1, length(vec))){
#    p = var_func(vec[i])
#    if (p >= pmin & p <= pmax){
#      bools[[i]] = T
#    }
#    else{
#      bools[[i]] = F
#    }
#  }
#  bools = unlist(bools)
#  return(bools)
#}
##flt_bools = filter_by_percentile(rnaseq_variance[, 1], 0.999, 1)
##flt_genes = rownames(rnaseq_variance)[flt_bools]
##flt_rnaseq = rnaseq[,flt_genes]

#flt_rnaseq = rnaseq[,order(rnaseq_variance, decreasing=T)[1:20]]

#hgnc_list <- getBM(filters="ensembl_gene_id",
#                   attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                   values=colnames(flt_rnaseq),
#                   uniqueRows = F,
#                   verbose = 104,
#                   mart=mart)

#
#print("renaming now")
#colnames(flt_rnaseq) = hgnc_list$hgnc_symbol

#

###### Plot
#feat = 'rnaseq_flt'
#feat_data = flt_rnaseq

## save the heatmaps to separate pdf pages
#fn = paste0(results_folder, 'rnaseq_flt.pdf')
#pdf(fn, onefile = TRUE)
#plot_list = make_heatmap(feat, feat_data, out,
#                         vis_fill_vals = T,
#                         clust_row=F, 
#                         clust_cols = F)
#dev.off()

#

#####################################################
## Plot Olink genes with most variance #
#####################################################

## filtering based on the variances of gene expression
## filtering is done using a percentile min and max
#olink_variance = olink %>% summarise_all(var)
#olink_variance =  t(olink_variance)

#flt_olink = olink[,order(olink_variance, decreasing=T)[1:20]]

#
##print("renaming now")
#olink_meta = read.table("../cmi-pb-pertussis/output/database_dump/olink_prot_info_dump.tsv",
#                              header=T)

#hgnc_list <- getBM(filters="uniprot_gn_id",
#                   attributes=c("uniprot_gn_id", "hgnc_symbol"),
#                   values=olink_meta$uniprot_id,
#                   uniqueRows = F,
#                   verbose = 104,
#                   mart=mart)

#olink_meta = merge(olink_meta, 
#                   hgnc_list,
#                   by.x = 'uniprot_id',
#                   by.y='uniprot_gn_id')

#olink_meta = olink_meta[!duplicated(olink_meta), ]

## rename the columns of the olink dataframe
#final_hgnc = match(colnames(flt_olink), olink_meta$olink_id)
#final_hgnc = olink_meta$hgnc_symbol[final_hgnc]
#colnames(flt_olink) = final_hgnc

## remove na's 
#na_cols = is.na(colnames(flt_olink))
#olink = olink[, !na_cols]

###### Plot
#feat = 'olink_flt'

#feat_data = flt_olink

## save the heatmaps to separate pdf pages
#fn = paste0(results_folder, 'olink_flt.pdf')
#pdf(fn, onefile = TRUE)
#plot_list = make_heatmap(feat, feat_data, out,
#                         vis_fill_vals = T,
#                         clust_row=F, 
#                         clust_cols = F)
#dev.off()
