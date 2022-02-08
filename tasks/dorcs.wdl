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
    }
    
    command {
        papermill /dorcs_jplot_notebook.ipynb ${output_filename} 
        -p rnaCountMatrix ${rna_matrix} 
        -p atacFragFile ${atac_fragments} 
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
    }

    output {
        File output_tsv = output_filename
        Array[File] plots = glob("/plots/*.pdf")
    }

    runtime {
        docker: docker_image
    }

    meta {
        author: "Siddarth Wekhande"
    }
}


workflow dorcs_workflow {
    input {
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
    }
        
    call find_dorcs { 
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
            topNGene = topNGene
    }
    
}
