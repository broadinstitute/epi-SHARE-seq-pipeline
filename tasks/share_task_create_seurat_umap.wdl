version 1.0

task create_seurat_umap {
    meta {
        version: 'v0.1'
        author: 'Kevin Dong (kdong@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: create seurat umap task'
    }

    input {
        #This tasks takes in an RNA matrix file and creates a UMAP using Seurat functions

        File rna_matrix = "rna.h5"
        String papermill = "TRUE"
        String genome = "hg38"
        
        String output_filename = "output.ipynb"
        String docker_image = ""
        Int mem_gb = 16
    }

    command {
        
        papermill $(which umap_seurat.ipynb) ${output_filename} \
        -p rna_matrix ${rna_matrix} \
        -p genome ${genome} \
        -p papermill ${papermill}
    }

    output {
        File notebook_output = output_filename
        File seurat_violin_plot = "plots/violinPlot.png"
        File seurat_mitochondria_qc_plot = "plots/mitochondria.png"
        File seurat_features_plot = "plots/features.png"
        File seurat_PCA_dim_loadings_plot = "plots/dimLoadings.png"
        File seurat_PCA_plot = "plots/pca.png"
        File seurat_heatmap_plot = "plots/heatmap.png"
        File seurat_jackstraw_plot = "plots/jackstraw.png"
        File seurat_elbow_plot = "plots/elbow.png"
        File seurat_umap_plot = "plots/umap.png"
        File plots_zip = "plots.zip"
    }

    runtime {
        cpu : 4
        memory : mem_gb+'G'
        docker : docker_image
    }
    
    parameter_meta {
        rna_matrix: {
            description: 'RNA matrix h5',
            help: 'RNA counts in matrix .h5 format',
            example: 'rna.h5'
        }
        
        papermill: {
            description: 'Boolean papermill flag',
            help: 'Flag to notebook run in papermill mode',
            example: 'TRUE'
        }
        
        genome: {
            description: 'Reference name',
            help: 'The name genome reference used to align.',
            examples: ['hg38', 'mm10', 'hg19', 'mm9'],
        }
        
        output_filename: {
            description: 'Output jupyter notebook name',
            help: 'The name assigned to output jupyter notebook',
            examples: 'output.ipynb',
        }
        
        docker_image: {
            description: 'Docker image.',
            help: 'Docker image for preprocessing step. 
            example: ['put link to gcr or dockerhub']
        }
    }
}