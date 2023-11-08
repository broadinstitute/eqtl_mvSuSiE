version 1.0

workflow run_mvSuSiE {
    input {
        Array[String] finemapped_qlts # ex. ips_D0.SuSiE_summary_with_cis_nom.parquet
        String combined_covariates # universal covariates to regress out (ips D0) e.g. with_X_ips_D0.​5PEERs.​combined_covariates.​txt
        Array[String] sample_names # ex. ips_D0 hep_D2 hep_D4 ...
        Array[File] expression_beds # filepaths for expression beds, in same order as sample names
        File inferred_cov_pcs # PCs from pipeline pt2 inferred covs (TODO: Add step for calculating this)
        String plink_file_prefix # ex. gs://landerlab-vcfs/StanleyCenter_CIRM_iPSC_WGS_callset_2021_01/maf01/WGS.filtered.plink
        String annotation_gtf # ex. gs://landerlab-20210915-ssong-macrophage-eqtls/resources/gencode.v26.GRCh38.genes.collapsed_only.gtf
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/eqtl_mvsusie:0.2'
    }

    call get_genes {
        input:
            finemapped_qlts=finemapped_qlts,
            docker_image=docker_image
    }

    call get_genotype_variant_tables {
        input:

    }

    scatter (gene in get_genes.gene_list){
        call run_mvSuSiE {
            input:
                gene=gene,
                combined_covariates=combined_covariates,
                sample_names=sample_names,
                expression_beds=expression_beds,
                inferred_cov_pcs=inferred_cov_pcs,
                plink_file_prefix=plink_file_prefix,
                annotation_gtf=annotation_gtf,
                docker_image=docker_image
        }
    }

    output {

    }
}

task get_genes {
    input {
        Array[String] finemapped_qlts
        String docker_image
    }
    command {
        set -ex
        (git clone https://github.com/broadinstitute/eqtl_mvSuSiE.git /app ; cd /app)
        python /get_genes.py -q ${sep=' ' finemapped_qlts}
    }
    output {
        Array[String] genes_list = read_lines("list_of_genes.txt")
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}

task run_mvSuSiE{
    input {
        String gene
        String combined_covariates
        Array[String] sample_names
        Array[File] expression_beds
        File inferred_cov_pcs
        String plink_file_prefix
        String annotation_gtf
        String docker_image
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eqtl_mvSuSiE.git /app ; cd /app)
        gsutil -m cp ${plink_file_prefix}.bed my_plink.bed
        gsutil -m cp ${plink_file_prefix}.bim my_plink.bim
        gsutil -m cp ${plink_file_prefix}.fam my_plink.fam
        gsutil -m cp ${annotation_gtf} annotation_gtf.gtf
        mkdir expression_beds
        gsutil -m cp ${sep=' ' expression_beds} expression_beds
        # TODO: figure out how to read in expression beds nicely to the python file (names of files?)
        python /get_tensorqtl_susie_map.py ${gene} ${inferred_cov_pcs} my_plink annotation_gtf.gtf ${combined_covariates} -s ${sep=' ' sample_names} -e ${sep expression_beds}


    }
}