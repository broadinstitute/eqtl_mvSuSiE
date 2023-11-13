import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("gene", type=str)
    parser.add_argument("-p", dest="phenotype_files", nargs='+', default=[],
                        help="Array of strings (filenames) of phenotype dfs")
    args = parser.parse_args()

    phenotypes = []
    for phenotype_fname in args.phenotype_files:
        phenotype_df = pd.read_table(phenotype_fname, index_col=0).T
        phenotypes.append(phenotype_df)

    large_phenotype_df = pd.concat(phenotypes, axis=1)
    large_phenotype_df.to_csv(f'{args.gene}_tensorqtl_regressed_phenotypes.csv', sep='\t')
