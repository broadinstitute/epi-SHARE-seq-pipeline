version 1.0

task find_dorcs {
    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: find DORCs task'
    }

    input {
        #This task takes in the RNA and ATAC files and finds the DORCs based on the cut-off criteria provided

        File rna_matrix
        File atac_fragments
        File peak_file

        String genome = "hg38"
        Int n_cores = 4
        String save_plots_to_dir = "TRUE"
        String output_filename = "output.ipynb"

        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        Int bootstraps = 100

        String docker_image = "polumechanos/dorc_task_find_dorcs"
    }

    command {
        gzip -dc ${atac_fragments} > tmp_fragments.bedpe

        papermill /dorcs_jplot_notebook.ipynb ${output_filename}
        -p rnaCountMatrix ${rna_matrix}
        -p atacFragFile tmp_fragments.bedpe
        -p peakFile ${peak_file}
        -p savePlotsToDir ${save_plots_to_dir}
        -p nCores ${n_cores}
        -p genome ${genome}
        -p minFeature_RNA ${minFeature_RNA}
        -p maxFeature_RNA ${maxFeature_RNA}
        -p percentMT_RNA ${percentMT_RNA}
        -p minCells_RNA ${minCells_RNA}
        -p dorcGeneCutOff ${dorcGeneCutOff}
        -p fripCutOff ${fripCutOff}
        -p corrPVal ${corrPVal}
        -p topNGene ${topNGene}
        -p windowPadSize ${windowPadSize}
        -p bootstraps ${bootstraps}
    }

    output {
        File output_tsv = output_filename
        Array[File] plots = glob("./plots/*.pdf")
    }

    runtime {
        docker: docker_image
    }
}


