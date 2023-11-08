version 1.0

workflow run_mvSuSiE {
    input {
        Array[File] finemapped_qlts # ex. ips_D0.SuSiE_summary_with_cis_nom.parquet
        File inferred_cov_pcs # PCs from pipeline pt2 inferred covs (TODO: Add step for calculating this)
        String plink_file_prefix # ex. gs://landerlab-vcfs/StanleyCenter_CIRM_iPSC_WGS_callset_2021_01/maf01/WGS.filtered.plink
        String annotation_gtf # ex. gs://landerlab-20210915-ssong-macrophage-eqtls/resources/gencode.v26.GRCh38.genes.collapsed_only.gtf
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/eqtl_mvsusie:0.2'
    }

    call get_genes {
        input:
            finemapped_qlts=finemapped_qlts
    }



    output {

    }
}

task get_genes {
    input {
        Array[File] finemapped_qlts
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