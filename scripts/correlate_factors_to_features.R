library("ggplot2")
library("reshape2")
library("stringr")
library("ggpubr")
library("gridExtra")
library("dplyr")
library("ComplexHeatmap")
library("circlize")
library('biomaRt')

# load the factors 
load("results/factors/factorizations.RData")

cytof_fn = 'data/cmi_pb/cytof.pma.day0.proc.tsv'
cytof = read.table(cytof_fn, header=T, row.names = 1)

olink_fn = 'data/cmi_pb/olink.pma.day0.proc.tsv'
olink = read.table(olink_fn, header=T, row.names = 1)

rnaseq_fn = 'data/cmi_pb/rnaseq.pma.day0.proc.tsv'
rnaseq = read.table(rnaseq_fn, header=T, row.names = 1)

features = list(cytof=cytof, olink=olink, rnaseq=rnaseq)

results_folder = "results/factors_v_features/"
dir.create(results_folder, showWarnings = F)

# loading the biomart
ensembl = useMart("ensembl", host="useast.ensembl.org", version = 104)
mart <- useDataset("hsapiens_gene_ensembl", ensembl)

##################################
##### Preprocessing the data #####
##################################

# change the metagene ensembl id to hgnc name
for (i in seq(1, length(out$factorizations))){
  genes = rownames(out$factorizations[[i]][[2]][[1]])
  
  print("obtaining hgnc name")
  hgnc_list <- getBM(filters="ensembl_gene_id",
                     attributes=c("ensembl_gene_id", "hgnc_symbol"),
                     values=genes,
                     uniqueRows = F,
                     verbose = 104,
                     mart=mart)
  print("renaming now")
  new_factored_data =  out$factorizations[[i]][[2]][[1]]
  new_factored_data = merge(new_factored_data, hgnc_list, by.x=0, by.y='ensembl_gene_id')
  
  print("running distinct")
  new_factored_data = distinct(new_factored_data, hgnc_symbol, .keep_all = T)
  rownames(new_factored_data) = new_factored_data$hgnc_symbol
  
  print("removing columns")
  new_factored_data = subset(new_factored_data, select = -c(Row.names, hgnc_symbol))
  
  print("saving to a new list in the out variable")
  out$factorizations[[i]][[3]] = list()
  out$factorizations[[i]][[3]][[1]] = new_factored_data
  print("renaming complete")
}

############################################
#### Plot Heatmap Using ComplexHeatmaps ####
############################################

# Function: correlates components versus feature data
make_heatmap <- function(feat,
                         feat_data,
                         factor_out,
                         clust_cols=T,
                         clust_rows=T, 
                         vis_fill_vals=F){
  
  plot_list = list()
  
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
  
  # correlate the tasks and factors
  for (i in seq(1,5)){
    
    # harmonize the data
    method = factor_out$method[[i]]
    factors = factor_out$factorizations[[i]][[1]]
    row.names(factors) = gsub('X', '', row.names(factors))
    shared_samples = intersect(row.names(factors), row.names(feat_data))
    final_factors = factors[shared_samples, ]
    final_features = feat_data[shared_samples, ]
    
    # correlate the factors and task
    factor_corrs = cor(final_factors, final_features)
    factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
    rownames(factor_corrs) = gsub("comp|SynVar", '\\1', rownames(factor_corrs))
    rownames(factor_corrs) = paste0('F', rownames(factor_corrs))
    color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    # make a heatmap of the correlations
    title = sprintf('%s factors versus %s', method, toupper(feat))
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
    print(p)
  }
  return(plot_list)
}

for (feat_name in names(features)){
  
  print(feat_name)
  
  # save the heatmaps to separate pdf pages
  fn = paste0(results_folder, feat_name, '2.pdf')
  pdf(fn, onefile = TRUE)
  
  # Set cytof configuration or RNA-seq/Olink configuration
  clust_cols = T
  vis_fill_vals = F
  if (feat_name == 'cytof'){
    clust_cols = T
    vis_fill_vals = T
  }
  
  # generate heatmaps in ggplot
  plot_list = make_heatmap(feat_name, features[[feat_name]], out,
                           clust_row = F, 
                           clust_cols = clust_cols, 
                           vis_fill_vals = vis_fill_vals)
  # close pdf
  dev.off()
  
}

# Find genes with large correlations in at least one factor
# correlate the tasks and factors
factor_out = out
feat_name = 'cytof'
feat_data = features[[feat_name]]
for (i in seq(1,5)){
  
  # harmonize the data
  method = factor_out$method[[i]]
  factors = factor_out$factorizations[[i]][[1]]
  row.names(factors) = gsub('X', '', row.names(factors))
  shared_samples = intersect(row.names(factors), row.names(feat_data))
  final_factors = factors[shared_samples, ]
  final_features = feat_data[shared_samples, ]
  
  # correlate the factors and task
  factor_corrs = cor(final_factors, final_features)
  factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
  rownames(factor_corrs) = gsub("comp|SynVar", '\\1', rownames(factor_corrs))
  rownames(factor_corrs) = paste0('F', rownames(factor_corrs))
  
  summary = apply(factor_corrs,
                  2,
                  FUN=function(x) return(sum(abs(x) > 0.75)))
  print(summary)

}







####################################################
# Plot RNA-seq genes with most variance#
####################################################

# filtering based on the variances of gene expression
# filtering is done using a percentile min and max
rnaseq_variance = rnaseq %>% summarise_all(var)
rnaseq_variance =  t(rnaseq_variance)

# Returns a vector for easy filtering
# filter based on the percentile using pmin and pmax
filter_by_percentile <- function(vec, pmin, pmax){
  var_func = ecdf(rnaseq_variance[, 1])
  bools = list()
  for (i in seq(1, length(vec))){
    p = var_func(vec[i])
    if (p >= pmin & p <= pmax){
      bools[[i]] = T
    }
    else{
      bools[[i]] = F
    }
  }
  bools = unlist(bools)
  return(bools)
}
#flt_bools = filter_by_percentile(rnaseq_variance[, 1], 0.999, 1)
#flt_genes = rownames(rnaseq_variance)[flt_bools]
#flt_rnaseq = rnaseq[,flt_genes]

flt_rnaseq = rnaseq[,order(rnaseq_variance, decreasing=T)[1:20]]

hgnc_list <- getBM(filters="ensembl_gene_id",
                   attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   values=colnames(flt_rnaseq),
                   uniqueRows = F,
                   verbose = 104,
                   mart=mart)


print("renaming now")
colnames(flt_rnaseq) = hgnc_list$hgnc_symbol



##### Plot
feat = 'rnaseq_flt'
feat_data = flt_rnaseq

# save the heatmaps to separate pdf pages
fn = paste0(results_folder, 'rnaseq_flt.pdf')
pdf(fn, onefile = TRUE)
plot_list = make_heatmap(feat, feat_data, out,
                         vis_fill_vals = T,
                         clust_row=F, 
                         clust_cols = F)
dev.off()



####################################################
# Plot Olink genes with most variance#
####################################################

# filtering based on the variances of gene expression
# filtering is done using a percentile min and max
olink_variance = olink %>% summarise_all(var)
olink_variance =  t(olink_variance)

flt_olink = olink[,order(olink_variance, decreasing=T)[1:20]]


#print("renaming now")
olink_meta = read.table("../cmi-pb-pertussis/output/database_dump/olink_prot_info_dump.tsv",
                              header=T)

hgnc_list <- getBM(filters="uniprot_gn_id",
                   attributes=c("uniprot_gn_id", "hgnc_symbol"),
                   values=olink_meta$uniprot_id,
                   uniqueRows = F,
                   verbose = 104,
                   mart=mart)

olink_meta = merge(olink_meta, 
                   hgnc_list,
                   by.x = 'uniprot_id',
                   by.y='uniprot_gn_id')

olink_meta = olink_meta[!duplicated(olink_meta), ]

# rename the columns of the olink dataframe
final_hgnc = match(colnames(flt_olink), olink_meta$olink_id)
final_hgnc = olink_meta$hgnc_symbol[final_hgnc]
colnames(flt_olink) = final_hgnc

# remove na's 
na_cols = is.na(colnames(flt_olink))
olink = olink[, !na_cols]

##### Plot
feat = 'olink_flt'

feat_data = flt_olink

# save the heatmaps to separate pdf pages
fn = paste0(results_folder, 'olink_flt.pdf')
pdf(fn, onefile = TRUE)
plot_list = make_heatmap(feat, feat_data, out,
                         vis_fill_vals = T,
                         clust_row=F, 
                         clust_cols = F)
dev.off()
