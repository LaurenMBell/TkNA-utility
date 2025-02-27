#
# loads in a csv file with gene pairs, and turns each into a scatterplot
# for each model in the data (LPS, DSS, and VECPAC)
#

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import argparse
import re
import csv

#sets the base directory for the data files
BASE_DIR = Path("~/TkNA/example_datasets_and_commands/tkna_input_deseq_normalized")

def load_expression_data(model):
    #loads in the gene expression data
    filepath = BASE_DIR / f"{model}.csv"
    try:
        data = pd.read_csv(filepath)
        return data
    except FileNotFoundError:
        print(f"Error: {filepath} not found :( ")
        return None

def load_group_map(model):
    #loads in map file
    group_path = BASE_DIR / f"{model}_group_map.csv"
    try:
        group_map = pd.read_csv(group_path, header=None, names=['ID', 'group'], sep=",")
        print(f"Loaded in the group map for {model}!")
        #print(group_map.head())
        return group_map
    except FileNotFoundError:
        print(f"Error: {group_path} not found :(")
        return None

def load_gene_pairs(filepath):
    pair_path = Path(filepath)
    
    # Check if the file exists
    if not pair_path.exists():
        print(f"Uh-oh: Gene pairs file not found at {pair_path} :(")
        return []
    
    gene_pairs = []
    
    #read in the gene pairs from the given csv file
    with open(pair_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                gene_pairs.append((row[0].strip(), row[1].strip()))
            elif len(row) == 1 and ',' in row[0]:
                genes = row[0].split(',')
                if len(genes) >= 2:
                    gene_pairs.append((genes[0].strip(), genes[1].strip()))
        
    print(f"Loaded gene pairs from {pair_path}!")
        
    return gene_pairs

def get_base_gene_id(gene_id):
    #extract the base ENSMUSG ID without any suffix
    match = re.match(r'(ENSMUSG\d+)', gene_id)
    return match.group(1) if match else gene_id

def create_scatterplot(gene1, gene2, data, groupmap, model, output_dir):
    #makes the scatterplot for each gene pair
    output_dir.mkdir(exist_ok=True, parents=True)
    
    control_label = 'treatment' if model == 'VECPAC' else 'control'
    title_group = 'CRE-' if model == 'VECPAC' else 'CTRL'
    
    control_samples = groupmap[groupmap['group'] == control_label]['ID'].astype(str).tolist()
    
    # extract base gene IDs for lookup
    base_gene1 = get_base_gene_id(gene1)
    base_gene2 = get_base_gene_id(gene2)
    
    # look for genes that match the base ID
    gene1_row = data[data.iloc[:, 0].str.contains(base_gene1)]
    gene2_row = data[data.iloc[:, 0].str.contains(base_gene2)]
    
    # If multiple matches found, take the first one
    if len(gene1_row) > 1:
        gene1_row = gene1_row.iloc[0:1]
    if len(gene2_row) > 1:
        gene2_row = gene2_row.iloc[0:1]
    
    if gene1_row.empty or gene2_row.empty:
        print(f"Bad news: couldn't find gene {gene1 if gene1_row.empty else gene2} in the data...")
        return None
    
    control_indices = [col for col in data.columns if str(col) in control_samples]
    
    gene1_values = gene1_row[control_indices].values.flatten()
    gene2_values = gene2_row[control_indices].values.flatten()
    
    plt.figure(figsize=(8, 6))
    plt.scatter(gene1_values, gene2_values)
    
    #I'll be honest, one of the graphs just WOULD NOT SHOW the regression line
    #no matter what I tried, so I used claude to debug and make this 
    # if statement (it worked tho, so I'll give it credit there)
    valid_indices = ~np.isnan(gene1_values) & ~np.isnan(gene2_values)
    if np.sum(valid_indices) >= 2:
        try:
            slope, intercept, r_value, _, _ = stats.linregress(gene1_values[valid_indices], gene2_values[valid_indices])

            x_min = min(gene1_values[valid_indices])
            x_max = max(gene1_values[valid_indices])
        
            x = np.linspace(x_min, x_max, 100)
            y = intercept + slope * x
        
            plt.plot(x, y, color='black', linewidth=1.5)
        
            plt.annotate(f'R² = {r_value**2:.3f}', xy=(0.05, 0.95), xycoords='axes fraction',
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
        except Exception as e:
            logger.error(f"Error calculating regression: {e}")
    
    plt.xlabel(f'{gene1}')
    plt.ylabel(f'{gene2}')
    plt.title(f'{gene1} vs. {gene2},\n{model} {title_group}')
    
    filename = f'{gene1}_vs_{gene2}_{model}_{title_group}.png'
    plt.savefig(output_dir / filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename

def process_gene_pairs(gene_pairs, output_dir):
    #makes the scatterplots for each model
    models = ['VECPAC', 'LPS', 'DSS']
    results = []
    
    model_data = {model: load_expression_data(model) for model in models}
    model_groupmaps = {model: load_group_map(model) for model in models}
    
    for gene1, gene2 in gene_pairs:
        pair_results = {'gene_pair': f'{gene1} vs {gene2}'}
        
        for model in models:
            if model_data[model] is not None and not model_data[model].empty and model_groupmaps[model] is not None and not model_groupmaps[model].empty:
                plot_filename = create_scatterplot(gene1, gene2, model_data[model], model_groupmaps[model], model, output_dir)
                pair_results[f'{model}_plot'] = plot_filename or 'Failed'
            else:
                pair_results[f'{model}_plot'] = 'Data not available'
        
        results.append(pair_results)
    
    return pd.DataFrame(results)

def main():
    #parse the arguments given and run the analysis and make the scatterplots
    parser = argparse.ArgumentParser()
    parser.add_argument('--gp', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()
    
    output_dir = Path(args.outdir).expanduser()
    output_dir.mkdir(exist_ok=True, parents=True)
    
    gene_pairs = load_gene_pairs(args.gp)
    
    results = process_gene_pairs(gene_pairs, output_dir)
    results.to_csv(output_dir / 'scatterplot_summary.csv', index=False)

    print(f"Scatterplots saved!")
    print("All done!")

if __name__ == "__main__":
    main()