version 1.0


# Import the tasks called by the pipeline
import "../tasks/dorcs_task_find_dorcs.wdl" as find_dorcs

workflow wf_dorcs {

    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org)'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to find DORCs from SHARE-seq data.'
    }

    input {
        File? rna_matrix
        File? atac_fragments
        File peak_file

        String genome
        Int n_cores = 4
        String save_plots_to_dir = "TRUE"
        String? output_filename

        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20
        Int windowPadSize = 50000

        Int numNearestNeighbor = 30
        Float numBackgroundPairs = 100000
        Float chunkSize = 50000

        String? prefix
        Int mem_gb = 64
        Int disk_gb = 100
        String? docker
    }
  
    File rna_matrix_ = select_first([rna_matrix])
    File atac_fragments_ = select_first([atac_fragments])

    if ( !defined(rna_matrix) || !defined(atac_fragments) ){
        call raise_exception as missing_input {
            input:
                  msg = "The genes-by-cell matrix or the dna fragments file are missing."
        }
    }

    call find_dorcs.find_dorcs as find_dorcs{
        input:
            rna_matrix = rna_matrix_,
            atac_fragments = atac_fragments_,
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
            docker_image = docker,
            prefix = prefix
    }

    output {
        File dorcs_notebook_output = find_dorcs.notebook_output
        File dorcs_notebook_log = find_dorcs.notebook_log
        File? seurat_violin_plot = find_dorcs.seurat_violin_plot
        File? j_plot = find_dorcs.j_plot
        File? plots_zip = find_dorcs.plots_zip
        File? dorcs_genes_summary = find_dorcs.dorcs_genes_summary
        File? dorcs_regions_summary = find_dorcs.dorcs_regions_summary
    }

}

# Task to report errors to user.
# From https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl
task raise_exception {
  input {
    String msg
    Array[String]? vals
  }
  command {
    echo -e "\n* Error: ${msg}\n" >&2
    echo -e "* Vals: ${sep=',' vals}\n" >&2
    exit 2
  }
  output {
    String error_msg = '${msg}'
  }
  runtime {
    maxRetries : 0
    cpu : 1
    memory : '2 GB'
    time : 1
    disks : 'local-disk 10 SSD'
  }
}
