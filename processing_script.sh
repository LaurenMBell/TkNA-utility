echo "process initiated!\n"

mkdir ./example_datasets_and_commands/choroid_data/output/histogram_output
echo "histogram output directory created!\n"

python ./reconstruction/intake_data.py --data-dir ./example_datasets_and_commands/choroid_data/choroid_plexus_data_deseq_normalized --out-file ./example_datasets_and_commands/choroid_data/output/all_data_and_metadata.zip
echo "data intake complete!\n"

# define the array of p-values
p_values=(0.05 0.1 0.15 0.2 0.25 0.3)

# loop through the p-values
for pval in "${p_values[@]}"; do
    echo "p-val threshold = $pval\n"

    # run the commands for each p-value
    python ./reconstruction/run.py --data-source ./example_datasets_and_commands/choroid_data/output/all_data_and_metadata.zip --config-file ./example_datasets_and_commands/choroid_data/choroid_plexus_data_deseq_normalized/config_corr_ind_pval_${pval}.json --out-file ./example_datasets_and_commands/choroid_data/output/network_output_pval_${pval}.zip
    python ./reconstruction/to_csv.py --data-file ./example_datasets_and_commands/choroid_data/output/network_output_pval_${pval}.zip --config-file ./example_datasets_and_commands/choroid_data/choroid_plexus_data_deseq_normalized/config_corr_ind_pval_${pval}.json --out-dir ./example_datasets_and_commands/choroid_data/output/network_output_pval_${pval}
    python ./visualization/corr_coef_histogram.py --data-file ./example_datasets_and_commands/choroid_data/output/network_output_pval_${pval}/network_output_comp.csv --pval $pval --out-dir ./example_datasets_and_commands/choroid_data/output/histogram_output

    echo "done with p-val $pval!\n\n"
done

# final step to combine all the histograms
python ./visualization/combine_histograms.py --input-dir ./example_datasets_and_commands/choroid_data/output/histogram_output --output-dir ./example_datasets_and_commands/choroid_data/output/histogram_output/final_figures
