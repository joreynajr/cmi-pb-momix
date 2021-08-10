library("ggplot2")
library("clusterCrit")
source("scripts/runfactorization.R")

# Performing the jDR's on the CMI-PB dataset
# The performances of the 5 jDR methods are here.

########## Data preprocessing ##########
# The data has been reformatted in Process_cmi_pb_data.R
# so as to be read by run_factorization function.

# Folder for results
results_folder <- "results/factors/"

# Create output folder
dir.create(results_folder, showWarnings = FALSE)

########## Running comparison ##########
# Two factor are then detected for each
# jDR method and the distribution of the cells
# with respect of Factor1 and Factor2 is plotted
# as a scatter plot. The obtained plots are available
# in the Results folder. The capability of the
# different jDR methods to cluster the cells according
# to their cell line of origin is finally evaluated
# through the C-index, whose value is reported in
# the Results folder.
# Run factorization methods
omic_fns <- c("rnaseq.momix.day0.input.tsv",
              "olink.momix.day0.input.tsv",
              "cytof.momix.day0.input.tsv")
out <- runfactorization("results/input_data/",
                        omic_fns,
                        10, 
                        sep=" ", 
                        filtering="stringent")

# Save the factors
fn = paste0(results_folder, "factorizations.RData")
save.image(file = fn)

## For each factorization method
# Parameters for the plots
# dot_size <- 1.5
# dot_alpha <- 1.0
# xlabel <- "Factor 1"
# ylabel <- "Factor 2"
#for(i in 1:length(out$factorizations)){
#    
#    # Get factorization result
#    factors <- out$factorizations[[i]][[1]]

#    # Delete NAs
#    factors <- factors[!is.na(factors[,1]) & !is.na(factors[,2]), ]
#    sample_annot <- sample_annot[!is.na(sample_annot[,1]) & !is.na(sample_annot[,2]), ]

#    # Data to be plotted
#    df <- data.frame(x =  factors[,1], y = factors[,2], color_by = sample_annot[,2])
#    
#    # Plot results
#    p <- ggplot(df, aes_string(x = "x", y = "y")) + 
#       geom_point(aes_string(color = "color_by"), size=dot_size, alpha=dot_alpha) + 
#       xlab(xlabel) + ylab(ylabel) +
#       # scale_shape_manual(values=c(19,1,2:18)[seq_along(unique(shape_by))]) +
#       theme(plot.margin = margin(20, 20, 10, 10), 
#             axis.text = element_text(size = rel(1), color = "black"), 
#             axis.title = element_text(size = 16), 
#             axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
#             axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
#             axis.line = element_line(color = "black", size = 0.5), 
#             axis.ticks = element_line(color = "black", size = 0.5),
#             panel.border = elemen_blank(), 
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(), 
#             panel.background = element_blank(),
#             legend.key = element_rect(fill = "white"),
#             legend.text = element_text(size = 16),
#             legend.title = element_text(size =16)
#       )
#    p + scale_color_manual(values=c("#0072B2", "#D55E00", "#CC79A7"))
#    
#    # Export plot as JPEG image
#    ggsave(paste0(results_folder, "plot_",out$method[i],".jpg"))
#}
