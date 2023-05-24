version 1.0

task find_dorcs {
    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: find DORCs task'
    }

    input {
        #This task takes in the RNA and ATAC files and finds the DORCs based on the cut-off criteria provided

        #DORCs parameters
        File rna_matrix
        File atac_fragments
        File? peak_file
        String genome
        Int n_cores = 4
        String save_plots_to_dir = "TRUE"
        String prefix = "prefix"

        #RNA QC parameters
        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        #ATAC QC parameter
        Float fripCutOff = 0.3
        Float chunkSize = 50000

        #Background correlation parameters
        Int numNearestNeighbor = 100
        Float numBackgroundPairs = 100000

        #DORC genes parameter
        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        Int dorcGeneCutOff = 10
        Float corrPVal = 0.05
        Int topNGene = 20

        String output_filename = "${prefix}.dorcs.notebook.${genome}.ipynb"
        String docker_image = "us.gcr.io/buenrostro-share-seq/dorcs_task_find_dorcs"
        #String docker_image = "swekhande/shareseq-prod:share-task-dorcs"
        Int mem_gb = 64
        Int disk_gb = 100
    }

    #Output filepaths

    String violin_plot = '${prefix}.dorcs.plots.${genome}/${prefix}.dorcs.rna_violin_plot.${genome}.png'
    String jplot = '${prefix}.dorcs.plots.${genome}/${prefix}.dorcs.jplot.${genome}.png'
    String dorc_genes_summ = '${prefix}.dorcs.dorc_genes_summary.${genome}.csv'
    String all_regions_summ = '${prefix}.dorcs.all_regions_summary.${genome}.csv'
    String plots_zip_dir = '${prefix}.dorcs.plots.${genome}.zip'
    #String papermill_log_filename = 'papermill.logfile.txt'
    String log_filename = "log/${prefix}.dorcs.logfile.${genome}.txt"

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
        -p chunkSize ${chunkSize} \
        -p prefix ${prefix}
    }

    output {
        File notebook_output = output_filename
        File notebook_log = log_filename
        #File papermill_log = papermill_log_filename

        File? seurat_violin_plot = violin_plot
        File? j_plot = jplot
        File? plots_zip = plots_zip_dir

        File? dorcs_genes_summary = dorc_genes_summ
        File? dorcs_regions_summary = all_regions_summ


    }

    runtime {
        cpu : 4
        memory : mem_gb+'G'
        docker : docker_image
        disks : 'local-disk ${disk_gb} LOCAL'
        maxRetries : 0
    }
}


