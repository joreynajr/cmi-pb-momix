echo hallmarks
Rscript scripts/rnaseq_pathway_analysis.R \
            --path-name hallmarks \
            --path-database data/bio_annotations/MSigDB_v7.4/h.all.v7.4.symbols.gmt \
            > logs/rnaseq_pathway_analysis.hallmarks.out \
            2> logs/rnaseq_pathway_analysis.hallmarks.err

echo "c2_curated"
Rscript scripts/rnaseq_pathway_analysis.R \
            --path-name c2_curated \
            --path-database data/bio_annotations/MSigDB_v7.4/c2.all.v7.4.symbols.gmt \
            > logs/rnaseq_pathway_analysis.c2_curated.out \
            2> logs/rnaseq_pathway_analysis.c2_curated.err
#
#skipping, don't produce any results
#echo "c3_reg_target"
#Rscript scripts/rnaseq_pathway_analysis.R \
#            --path-name c3_reg_target \
#            --path-database data/bio_annotations/MSigDB_v7.4/c3.all.v7.4.symbols.gmt \
#            > logs/rnaseq_pathway_analysis.c3_reg_target.out \
#            2> logs/rnaseq_pathway_analysis.c3_reg_target.err
#
echo "c5_ontologies"
Rscript scripts/rnaseq_pathway_analysis.R \
            --path-name c5_ontologies \
            --path-database data/bio_annotations/MSigDB_v7.4/c5.all.v7.4.symbols.gmt \
            > logs/rnaseq_pathway_analysis.c5_ontologies.out \
            2> logs/rnaseq_pathway_analysis.c5_ontologies.err

#echo "c7_immunology"
#Rscript scripts/rnaseq_pathway_analysis.R \
#            --path-name c7_immunology \
#            --path-database data/bio_annotations/MSigDB_v7.4/c7.all.v7.4.symbols.gmt \
#            > logs/rnaseq_pathway_analysis.c7_immunology.out \
#            2> logs/rnaseq_pathway_analysis.c7_immunology.err

echo "c8_cell_signatures"
Rscript scripts/rnaseq_pathway_analysis.R \
            --path-name c8_cell_signatures \
            --path-database data/bio_annotations/MSigDB_v7.4/c8.all.v7.4.symbols.gmt \
            > logs/rnaseq_pathway_analysis.c8_cell_signatures.out \
            2> logs/rnaseq_pathway_analysis.c8_cell_signatures.err
