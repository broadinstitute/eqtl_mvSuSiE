import argparse
import pandas as pd
from collections import Counter

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", dest="qtl_finemaps", nargs='+', default=[],
                        help="Array of strings (filenames) of qtl finemapping results")
    args = parser.parse_args()

    high_pip_genes = list()
    for qtl_finemap_fname in args.qtl_finemaps:
        qtl_finemap = pd.read_parquet(qtl_finemap_fname)
        high_pip_finemap = qtl_finemap.query('pip > 0.5')
        high_pip_genes.extend(high_pip_finemap.phenotype_id.unique())

    lst = Counter(high_pip_genes)
    genes = [gene for gene, count in lst.items() if (count >= 3 and count < 6)]

    with open('list_of_genes.txt', mode='wt', encoding='utf-8') as fh:
        # source https://stackoverflow.com/questions/10188453/python-insert-line-break-every-x-number-of-elements
        for i, m in enumerate(genes, 1):
            fh.write(m + [' ', '\n'][i % 5 == 0])
