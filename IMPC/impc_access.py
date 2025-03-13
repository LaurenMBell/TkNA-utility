import csv
import requests
import argparse
from pathlib import Path


def fetch_gene_data(gene_symbol):
    #gets all the data from the IMPC database for a gene_symbol
    base_url = "https://www.ebi.ac.uk/mi/impc/solr/experiment/select"
    params = {
        "q": f"gene_symbol:{gene_symbol}",
        "wt": "json", #the API returns in json data
        "rows": 1000 #i could add more if needed 
    }
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        return response.json().get("response", {}).get("docs", [])
    else:
        print("THere was an error getting the data")
        return None


def data_to_csv(data, output_dir, gene_symbol):
    #saves the gene data to a csv file 
    gene_dir = Path(output_dir) / gene_symbol
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = gene_dir / f"{gene_symbol}_data.csv"
    
    fields = set()
    for item in data:
        fields.update(item.keys())
    
    with open(output_file, "w", newline="") as outfile:
        w = csv.DictWriter(outfile, fieldnames=sorted(fields))
        w.writeheader()
        w.writerows(data)
    
    print(f"Data written to {output_file}")
    return output_file


def process_gene_list(gene_list, output_dir):
    #processes all the genes in the given csv list
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    #read the list in
    genes = []
    with open(gene_list, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:
                genes.append(row[0].strip())

    
    #process each gene
    results = []
    for gene in genes:
        data = fetch_gene_data(gene)
        if data:
            output_file = data_to_csv(data, output_dir, gene)
            results.append({"gene": gene, "file": str(output_file), "records": len(data)})
        else:
            print(f"No data found for {gene}")
            results.append({"gene": gene, "file": "N/A", "records": 0})
    
    #write it all to a file
    return_file = output_path / "output.csv"
    with open(return_file, "w", newline="") as outfile:
        w = csv.DictWriter(outfile, fieldnames=["gene", "file", "records"])
        w.writeheader()
        w.writerows(results)
    
    print(f"All done!")


def main():
    #parse the arguments and run the processing 
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_list")
    parser.add_argument("--output_dir")
    args = parser.parse_args()
    
    process_gene_list(args.gene_list, args.output_dir)


if __name__ == "__main__":
    main()