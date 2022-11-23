library(ggplot2)

plot_list = list()
gene_sets = c("hallmarks", "c2_curated", #"c3_reg_target",
              "c5_ontologies", "c7_immunology", "c8_cell_signatures")

i = 1
for (gene_set in gene_sets){
    fn = sprintf("results/pathway_analysis/%s/report.tsv", gene_set)
    data = read.table(fn)
    data[is.na(data)] = 0
    data['method'] = rownames(data)

    # add jitter
    set.seed(0)
    data$nonZeroFacs = data$nonZeroFacs +  rnorm(nrow(data), 0, sd=0.075)
    
    print(gene_set)
    print(data)
    
    title = sprintf("Biological enrichment: %s", gene_set)
    title = gsub("_", " ", title)
    xlabel = "Number of metagenes enriched in at\nleast a signature from %s"
    xlabel = sprintf(xlabel, gene_set)
    p = ggplot(data, aes(x=nonZeroFacs, y=selectivity)) + 
            geom_point(aes(color=method), size=5, shape=18) + 
            ggtitle(title) + 
            scale_x_continuous(xlabel, breaks=seq(0,10,1), limits=c(-0.5,10.5)) +
            scale_y_continuous(breaks=seq(0.5,1,0.1), limits=c(0.45,1.025))
    
    plot_list[[i]] = p
    i = i + 1 
}

# save the heatmaps to separate pdf pages
fn = paste0("results/pathway_analysis/metagenes_v_selectivity.pdf")
theme_set(theme_grey(base_size=6))
pdf(fn, onefile = TRUE, width=3, height = 3)
for (i in seq(1, length(plot_list))){
    print(plot_list[[i]])
}
dev.off()
