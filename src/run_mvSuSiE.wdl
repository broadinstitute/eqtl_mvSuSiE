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
        File mashr_strong_prior
        String docker_image_py = 'us.gcr.io/landerlab-atacseq-200218/eqtl_mvsusie:0.5'
        String docker_image_r = 'us.gcr.io/landerlab-atacseq-200218/eqtl_mvsusie_r:0.2'
    }

    call get_genes {
        input:
            finemapped_qlts=finemapped_qlts,
            docker_image=docker_image_py
    }

    scatter (gene in get_genes.genes_list){
        call run_qtl_susie_regression {
            input:
                gene=gene,
                combined_covariates=combined_covariates,
                sample_names=sample_names,
                expression_beds=expression_beds,
                inferred_cov_pcs=inferred_cov_pcs,
                plink_file_prefix=plink_file_prefix,
                annotation_gtf=annotation_gtf,
                docker_image=docker_image_py
        }

        call concat_phenotypes {
            input:
                phenotype_files=run_qtl_susie_regression.phenotype_files,
                gene=gene,
                docker_image=docker_image_py
        }

        call run_mvSuSiE {
            input:
                gene=gene,
                mashr_strong_prior=mashr_strong_prior,
                genotype_file=run_qtl_susie_regression.genotype_file,
                phenotype_file=concat_phenotypes.phenotype_file,
                sample_names=sample_names,
                docker_image=docker_image_r
        }
    }

    call gather_mvsusie_outputs {
        input:
            mvsusie_results=run_mvSuSiE.mvsusie_results,
            docker_image=docker_image_py
    }

    output {
        Array[String] gene_list = get_genes.genes_list
        Array[File] mvsusie_results = run_mvSuSiE.mvsusie_results
        File mvsusie_results_all = gather_mvsusie_outputs.mvsusie_results_all
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
        python /app/src/get_genes.py -q ${sep=' ' finemapped_qlts}
        head -5 list_of_genes.txt > only_five_genes.txt
    }
    output {
        Array[String] genes_list = read_lines("only_five_genes.txt")
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}

task run_qtl_susie_regression {
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
        mkdir expression_beds_dir
        gsutil -m cp ~{sep=" " expression_beds} expression_beds_dir
        basenames=()
        for f in ~{sep=" " expression_beds}
        do
            basenames+=,expression_beds_dir/$(basename $f)
        done
        python /app/src/get_tensorqtl_susie_map.py ${gene} ${inferred_cov_pcs} my_plink annotation_gtf.gtf ${combined_covariates} -s ${sep=' ' sample_names} -e $basenames

    }

    output {
        Array[File] phenotype_files = glob("${gene}_tensorqtl_regressed_exp_*")
        File genotype_file = "${gene}_tensorqtl_regressed_genotypes.csv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}

task concat_phenotypes {
    input {
        Array[File] phenotype_files
        String gene
        String docker_image
    }
    command {
        set -ex
        (git clone https://github.com/broadinstitute/eqtl_mvSuSiE.git /app ; cd /app)
        python /app/src/concat_phenotypes.py ${gene} -p ${sep=' ' phenotype_files}
    }
    output {
        File phenotype_file = "${gene}_tensorqtl_regressed_phenotypes.csv"
    }

    runtime {
            docker: docker_image
            cpu: 1
            memory: "16GB"
            preemptible: 1
    }
}

task run_mvSuSiE {
    input {
        String gene
        File mashr_strong_prior
        File genotype_file
        File phenotype_file
        Array[String] sample_names
        String docker_image
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/eqtl_mvSuSiE.git /app ; cd /app)
        Rscript /app/src/run_mvSuSiE.R ${gene} ${mashr_strong_prior} ${genotype_file} ${phenotype_file} ${sep=',' sample_names}
        tr ',' '\t' < ${gene}_mvsusie_final_output.csv > ${gene}_mvsusie_final_output.tsv
    }
    output {
        File mvsusie_results = "${gene}_mvsusie_final_output.tsv"
    }
    runtime {
            docker: docker_image
            cpu: 1
            memory: "16GB"
            preemptible: 1
    }
}

task gather_mvsusie_outputs {
    input {
        Array[String] mvsusie_results
        String docker_image
    }
    String first_file = basename(mvsusie_results[0])

    command {
        set -ex
        mkdir mvsusie_results_dir
        gsutil -m cp ~{sep=" " mvsusie_results} mvsusie_results_dir
        head -n 1 mvsusie_results_dir/${first_file} > mvsusie_final_output_all.tsv
        tail -n +2 -q mvsusie_results_dir/* >> mvsusie_final_output_all.tsv
    }
    output {
        File mvsusie_results_all = "mvsusie_final_output_all.tsv"
    }
    runtime {
            docker: docker_image
            cpu: 1
            memory: "16GB"
            preemptible: 1
    }
}