rscript="/home/jreyna/software/anaconda3/envs/Momix_R-4.1.1/bin/Rscript"
PATH="/home/jreyna/software/anaconda3/envs/Momix_R-4.1.1/bin/":$PATH
# Runner.sh has been make to run each individual script 
# in topological ordering (DAG order)

# preprocess the data for momix scripts
read -r cmd << EOM
$rscript scripts/process_data.R \
               > logs/process_data.out \
               2> logs/process_data.err
EOM
#echo "Running $cmd"
#eval $cmd
#vim -p logs/process_data.out logs/process_data.err


# calculate factors
read -r cmd << EOM
$rscript scripts/calculate_factors.R \
               > logs/calculate_factors.out \
               2> logs/calculate_factors.err
EOM
echo "Running $cmd"
eval $cmd
vim -p logs/calculate_factors.out logs/calculate_factors.err


## Correlate Factors To Taskss
#read -r cmd << EOM
#$rscript scripts/correlate_factors_to_tasks.R \
#               > logs/correlate_factors_to_tasks.out \
#               2> logs/correlate_factors_to_tasks.err
#EOM
#
#echo "Running $cmd"
#eval $cmd
#vim -p logs/correlate_factors_to_tasks.out logs/correlate_factors_to_tasks.err

#$rscript scripts/correlate_factors_to_features.R \
#               > logs/correlate_factors_to_features.out \
#               2> logs/correlate_factors_to_features.err
#$rscript scripts/rnaseq_pathway_analysis.R > logs/rnaseq_pathway_analysis.out 2> logs/rnaseq_pathway_analysis.err
