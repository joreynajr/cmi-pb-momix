# Runner.sh has been make to run each individual script 
# in topological ordering (DAG order)
#Rscript scripts/process_data.R \
#               > logs/process_data.out \
#               2> logs/process_data.err
#Rscript scripts/calculate_factors.R \
#               > logs/calculate_factors.out \
#               2> logs/calculate_factors.err
#Rscript scripts/correlate_factors_to_tasks.R \
#               > logs/correlate_factors_to_tasks.out \
#               2> logs/correlate_factors_to_tasks.err
#Rscript scripts/correlate_factors_to_features.R \
#               > logs/correlate_factors_to_features.out \
#               2> logs/correlate_factors_to_features.err
Rscript scripts/rnaseq_pathway_analysis.R \
                > logs/rnaseq_pathway_analysis.out \
                2> logs/rnaseq_pathway_analysis.err