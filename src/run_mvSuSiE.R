args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
    library(mashr)
    library(mvsusieR)
    library(tidyverse)
    library(dplyr)
    library(susieR)
    library(Rfast)
})

message("Reading in arguments.")

# the mash object fit to the strong data, e.g. mashr_strong_mash_object.rds
mash_prior_positives <- readRDS(args[1])
# group_names (ips_d0, etc)
group_names <- unlist(strsplit(args[2], ","))
# genes
gene_names <- args[3:length(args)]

for (gene_name in gene_names){
    message(gene_name)
    # regressed genotypes from python generation,
    #   e.g. "HES4_tensorqtl_regressed_genotypes.csv"
    genotypes <- paste0("genotype_files_dir/", gene_name,
        "_tensorqtl_regressed_genotypes.csv")
    # regressed phenotypes, e.g. "HES4_tensorqtl_regressed_phenotypes.csv"
    phenotypes <- paste0("phenotype_files_dir/", gene_name,
        "_tensorqtl_regressed_phenotypes.csv")

    # read in the one regressed genotype matrix
    x <- read.table(genotypes, header = TRUE,
            row.names = 1, check.names = FALSE)
    y_mvsusie <- read.table(phenotypes, header = TRUE,
            row.names = 1, check.names = FALSE)

    x_mvsusie <- t(x[, rownames(y_mvsusie)])

    # generate the prior from the mash strong object
    message("Generating prior objects.")
    newUlist <- mashr:::expand_cov(mash_prior_positives$fitted_g$Ulist,
                                mash_prior_positives$fitted_g$grid,
                                usepointmass = TRUE)
    prior_mixt <- create_mixture_prior(list(matrices = newUlist[-1],
                                weights = mash_prior_positives$fitted_g$pi[-1]),
                                null_weight = mash_prior_positives$fitted_g$pi[1])


    message("Running mvSuSiE.")
    # run the mvsusie model
    fit_mv <- mvsusie(x_mvsusie, y_mvsusie, standardize = TRUE,
            prior_variance = prior_mixt,
            estimate_prior_variance = TRUE)

    message("Generating dataframe.")
    mvsusie_data_frame <- data.frame ()
    if (is.null(fit_mv$sets$cs)){
        data <- c('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')
        mvsusie_data_frame <- rbind(mvsusie_data_frame,
            c(gene_name, 'mvSusie', data))
    } else {
        for (j in 1:length(fit_mv$sets$cs)){
            L = fit_mv$sets$cs[[j]]
            credible_set <- fit_mv$sets$cs_index[j]
            for (variant_id in L){
                data <- c(colnames(x_mvsusie)[as.numeric(variant_id)],
                                fit_mv$pip[colnames(x_mvsusie)[as.numeric(variant_id)]],
                                credible_set)
                mvsusie_data_frame <- rbind(mvsusie_data_frame,
                                    c(gene_name, 'mvSusie', data,
                                    fit_mv$b1[credible_set, variant_id, ]))
            }
        }
    }
    colnames(mvsusie_data_frame) <- c('phenotype_id', 'Method', 'Variant_id', 'PIP', 'Credible Set', group_names)

    write.csv(mvsusie_data_frame, paste0(gene_name, "_mvsusie_final_output.csv"),
        row.names=FALSE)
    message("Done.")
}
