import argparse
import pandas as pd
import numpy as np
import torch
import tensorqtl
import qtl
from qtl import annotation
from tensorqtl import genotypeio, susie, core


def my_map_susie(
    variant_df,
    genotype_df,
    phenotype_df,
    phenotype_pos_df,
    covariates_df,
    window=1000000,
    maf_threshold=0.05,
):
    # reduced function from tensorql
    # TODO: get link for source
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    residualizer = core.Residualizer(
        torch.tensor(covariates_df.values, dtype=torch.float32).to(device)
    )
    # print(phenotype_df)
    genotype_ix = np.array(
        [genotype_df.columns.tolist().index(i) for i in phenotype_df.columns]
    )
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)

    igc = genotypeio.InputGeneratorCis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df, window=window
    )
    if igc.n_phenotypes == 0:
        raise ValueError("No valid phenotypes found.")

    copy_keys = ["pip", "sets", "converged", "elbo", "niter", "lbf_variable"]
    susie_summary = []

    for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(
        igc.generate_data(verbose=True), 1
    ):
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:, genotype_ix_t]
        # print(genotypes_t.shape)
        core.impute_mean(genotypes_t)

        variant_ids = variant_df.index[
            genotype_range[0] : genotype_range[-1] + 1
        ].rename("variant_id")

        # filter monomorphic variants
        mask_t = ~(genotypes_t == genotypes_t[:, [0]]).all(1)

        if maf_threshold > 0:
            maf_t = core.calculate_maf(genotypes_t)
            mask_t &= maf_t >= maf_threshold
        if mask_t.any():
            genotypes_t = genotypes_t[mask_t]
            mask = mask_t.cpu().numpy().astype(bool)
            variant_ids = variant_ids[mask]
            genotype_range = genotype_range[mask]

        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
        genotypes_res_t = residualizer.transform(genotypes_t)  # variants x samples
        phenotype_res_t = residualizer.transform(
            phenotype_t.reshape(1, -1)
        )  # phenotypes x samples

        ## RUN SUSIE
        res = susie.susie(
            genotypes_res_t.T,
            phenotype_res_t.T,
            L=10,
            scaled_prior_variance=0.2,
            coverage=0.95,
            min_abs_corr=0.5,
            estimate_residual_variance=True,
            estimate_prior_variance=True,
            tol=1e-3,
            max_iter=100,
        )

        af_t = genotypes_t.sum(1) / (2 * genotypes_t.shape[1])
        res["pip"] = pd.DataFrame(
            {"pip": res["pip"], "af": af_t.cpu().numpy()}, index=variant_ids
        )

        return res["pip"], phenotype_res_t, genotypes_res_t, genotype_range


def call_map_susie(
    gene_str,
    variant_df,
    genotype_df,
    universal_covariates,
    annot,
    group_names,
    expression_beds,
):
    """For a given gene:
        Get covariates, phenotypes, genotypes for every day.
        Run susie and get regressed results for mvsusie.
    Function akes a gene string, for example 'ABHD4', and the variant and genotype dfs generated in main.
    """
    phenotype_dfs = {}
    phenotype_pos_dfs = {}
    susie_res_dfs_mymap = {}

    phenotype_regr_dfs = {}
    genotype_regr_dfs = {}
    genotype_ranges = {}

    # the annotation give the location of the genes, and then you can query the genotype df
    gene = annot.get_gene(gene_str)  # pass in gene name
    # get all the variants within a window of a gene
    w = 1000000
    vt_df = variant_df[
        (variant_df["chrom"] == gene.chr)
        & (variant_df["pos"] >= gene.tss - w)
        & (variant_df["pos"] <= gene.tss + w)
    ]

    # genotype df is now universal across all the days, only includes the 181 donors that are present in each day
    gt_df = genotype_df.loc[vt_df.index, universal_covariates.index]

    for group_name, expression_bed in zip(group_names, expression_beds):
        # covariates now come from universal_covariates. this has the 5 PCs used across all days, sex age etc. info same across all days,
        # and the new 11 PCs taken from the PCA run on the inferred covariates of each day

        # expression matrices
        (
            phenotype_dfs[group_name],
            phenotype_pos_dfs[group_name],
        ) = tensorqtl.read_phenotype_bed(expression_bed)

        # filter the phenotype dfs ahead of time too to only be the donors that have data in every day
        phenotype_dfs[group_name] = phenotype_dfs[group_name][gt_df.columns]

        # cis permutation results
        # do we need this here? or just for picking genes
        # cis_dfs[group_name] = pd.read_parquet(f'gs://landerlab-20220124-ssong-village-eqtls/analysis_freeze/{group_name}/{group_name}.5PEERs.cis_qtl.sigificant.parquet')

        # fine-map
        (
            susie_res_dfs_mymap[group_name],
            phenotype_regr_dfs[group_name],
            genotype_regr_dfs[group_name],
            genotype_ranges[group_name],
        ) = my_map_susie(
            vt_df,
            gt_df,
            phenotype_dfs[group_name].loc[gene_str].to_frame().T,
            phenotype_pos_dfs[group_name].loc[gene_str].to_frame().T,
            universal_covariates,
            maf_threshold=0.05,
        )

        # Change from unlabeled numpy to pandas dfs with col & index labels
        phenotype_regr_dfs[group_name] = pd.DataFrame(
            phenotype_regr_dfs[group_name], columns=phenotype_dfs[group_name].columns
        )
        genotype_regr_dfs[group_name] = pd.DataFrame(
            genotype_regr_dfs[group_name],
            index=gt_df.iloc[genotype_ranges[group_name],].index,
            columns=phenotype_dfs[group_name].columns,
        )

    return susie_res_dfs_mymap, phenotype_regr_dfs, genotype_regr_dfs


def main():
    parser = argparse.ArgumentParser(
        description="For a given gene, get covariates, phenotypes, genotypes for every day. Run susie and get regressed results for mvsusie"
    )
    parser.add_argument("gene", type=str)
    parser.add_argument(
        "inferred_cov_pcs",
        type=str,
        help="File path for the PCs from the inferred covs",
    )
    parser.add_argument(
        "plink_prefix_path",
        type=str,
        help="General name for plink files, e.g. WGS.filtered.plink",
    )
    parser.add_argument(
        "annotation_gtf", type=argparse.FileType("r"), help="Gene chr gtf"
    )
    parser.add_argument(
        dest="combined_covariates",
        type=str,
        help="Combined covariates from pipeline pt2 to regress out",
    )
    parser.add_argument(
        "-s", dest="sample_names", nargs="+", default=[], help="Array of sample names"
    )
    parser.add_argument(
        "-e",
        dest="expression_beds",
        type=str,
        help="String of ',' separated expression bed files. Must be in same order as sample names",
    )
    args = parser.parse_args()

    # Get universal covariates to regress out
    print("Reading in covariates.")
    inferred_covs_pcs = pd.read_csv(args.inferred_cov_pcs, index_col=0)
    covariates_df = pd.read_csv(args.combined_covariates, sep="\t", index_col=0).T
    # TODO: Generalize for greater audience. for now we know we want ips_D0
    stagnant_covs = covariates_df.drop(
        covariates_df.filter(regex="Inferred").columns, axis=1
    )
    universal_covariates = inferred_covs_pcs.merge(
        stagnant_covs, left_index=True, right_index=True
    )

    # Get gene annotation names
    print("Getting gene annotation.")
    annot = qtl.annotation.Annotation(args.annotation_gtf.name)

    # Get genotype and variant dfs
    print("Getting genotype and phenotype info.")
    pr = genotypeio.PlinkReader(args.plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.index = (
        pr.bim["chrom"]
        + "_"
        + pr.bim["pos"].astype(str)
        + "_"
        + pr.bim["a1"]
        + "_"
        + pr.bim["a0"]
    )
    genotype_df.index = (
        pr.bim["chrom"]
        + "_"
        + pr.bim["pos"].astype(str)
        + "_"
        + pr.bim["a1"]
        + "_"
        + pr.bim["a0"]
    )

    # Delete duplicates
    genotype_df = genotype_df[~variant_df.duplicated(keep=False)]
    variant_df = variant_df[~variant_df.duplicated(keep=False)]

    print("Running susie & regressing.")
    # expression_beds files are read in as a ',' separated str (for ease because wdl is weird about bash variables)
    susie_res_dfs_mymap, phenotype_regr_dfs, genotype_regr_dfs = call_map_susie(
        args.gene,
        variant_df,
        genotype_df,
        universal_covariates,
        annot,
        args.sample_names,
        args.expression_beds.split(',')[1:],
    )

    print("Saving files.")
    for group_name in args.sample_names:
        phenotype_regr_dfs[group_name] = phenotype_regr_dfs[group_name].rename(index={0:group_name})
        phenotype_regr_dfs[group_name].to_csv(f'{args.gene}_tensorqtl_regressed_exp_{group_name}.csv', sep='\t')

    # all of the genotype regr dfs are the same now for every day
    genotype_regr_dfs['ips_D0'].to_csv(f'{args.gene}_tensorqtl_regressed_genotypes.csv', sep='\t')
    print("Done.")


if __name__ == "__main__":
    main()
