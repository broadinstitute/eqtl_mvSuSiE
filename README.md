# eqtl_mvSuSiE
Terra workflow for running mvSuSiE in R based on eQTL finemapping results from https://github.com/broadinstitute/eqtl_pipeline_terra.

## Description:
Takes in the prior from mash:
* from mashr, get canonical from strong matrices, and fit mash object
* Prior: Fit mash object to strong data & its canonical matrices. Call cov_expand to get the weights of every covariate matrix. Pass the weights and cov matrices to create_mixture_prior.

Takes in the PCA on inferred covariates.

Takes in finemapping & expression outputs from eqtl pipeline

Only runs mvsusie on genes that finemapped in >= 3 days but < 6 days

Runs mvsusie for 5 genes at a time, across all genes in parallel

Returns dataframe with columns 'phenotype_id', 'Method', 'Variant_id', 'PIP', 'Credible Set', and betas for each response variable (in our case, ipsd0, hepd2, etc).
Genes and the variants in each of their credible sets should have results.
