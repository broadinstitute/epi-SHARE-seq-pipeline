version 1.0


# Import the tasks called by the pipeline
import "../tasks/dorcs_task_find_dorcs.wdl" as find_dorcs

workflow wf_dorcs {
    input {
        File rna_matrix
        File atac_fragments
        File peak_file

        String genome
        Int n_cores = 4
        String save_plots_to_dir = "TRUE"
        String output_filename 

        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20
        Int windowPadSize = 50000
        
        Int numNearestNeighbor = 100
        Float numBackgroundPairs = 100000
        Float chunkSize = 50000
        
        Int mem_gb = 16
        Int disk_gb = 100
        String docker = "swekhande/shareseq-prod:share-task-dorcs"
    }

    call find_dorcs.find_dorcs as find_dorcs{
        input:
            rna_matrix = rna_matrix,
            atac_fragments = atac_fragments,
            peak_file = peak_file,
            genome = genome,
            n_cores = n_cores,
            save_plots_to_dir = save_plots_to_dir,
            output_filename = output_filename,
            minFeature_RNA = minFeature_RNA,
            maxFeature_RNA = maxFeature_RNA,
            percentMT_RNA = percentMT_RNA,
            minCells_RNA = minCells_RNA,
            dorcGeneCutOff = dorcGeneCutOff,
            fripCutOff = fripCutOff,
            corrPVal = corrPVal,
            topNGene = topNGene,
            windowPadSize = windowPadSize,
            numNearestNeighbor = numNearestNeighbor,
            numBackgroundPairs = numBackgroundPairs,
            chunkSize = chunkSize,
            mem_gb = mem_gb,
            disk_gb = disk_gb,
            docker_image = docker
    }

    output {
        File notebook_output = find_dorcs.notebook_output
        File? seurat_violin_plot = find_dorcs.seurat_violin_plot
        File? j_plot = find_dorcs.j_plot
        File? plots_zip = find_dorcs.plots_zip
        File? dorcs_genes_summary = find_dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = find_dorcs.dorcs_regions_summary
    }

}
