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

        String genome
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
        
        Int numNearestNeighbor = 100
        Float numBackgroundPairs = 100000
        Float chunkSize = 50000

        String docker_image = "swekhande/dorcs:dorcs-notebook"
        Int mem_gb = 16
        Int disk_gb = 100
    }

    command {
        gzip -dc ${atac_fragments} > tmp_fragments.bedpe

        papermill $(which dorcs_jplot_notebook.ipynb) ${output_filename} \
        -p rnaCountMatrix ${rna_matrix} \
        -p atacFragFile tmp_fragments.bedpe \
        -p peakFile ${peak_file} \
        -p savePlotsToDir ${save_plots_to_dir} \
        -p nCores ${n_cores} \
        -p genome ${genome} \
        -p minFeature_RNA ${minFeature_RNA} \
        -p maxFeature_RNA ${maxFeature_RNA} \
        -p percentMT_RNA ${percentMT_RNA} \
        -p minCells_RNA ${minCells_RNA} \
        -p dorcGeneCutOff ${dorcGeneCutOff} \
        -p fripCutOff ${fripCutOff} \
        -p corrPVal ${corrPVal} \
        -p topNGene ${topNGene} \
        -p windowPadSize ${windowPadSize} \
        -p numNearestNeighbor ${numNearestNeighbor} \
        -p numBackgroundPairs ${numBackgroundPairs} \
        -p chunkSize ${chunkSize}
    }

    output {
        File notebook_output = output_filename
        File seurat_violin_plot = "plots/RNAViolinPlot.png"
        File j_plot = "plots/JPlot.png"
        File plots_zip = "plots.zip"
        File dorcs_genes_summary = "dorc_genes_summary.csv"
        File dorcs_regions_summary = "dorc_regions_summary.csv"
    }

    runtime {
        cpu : 4
        memory : mem_gb+'G'
        docker : docker_image
        disks : 'local-disk ${disk_gb} LOCAL'
        maxRetries : 0
    }
}


