version 1.0

task create_archr_umap {
    meta {
        version: 'v0.1'
        author: 'Kevin Dong (kdong@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: create archr umap task'
    }

    input {
        #This tasks takes in an ATAC fragment file and creates a UMAP using ArchR functions

        File atac_frag
        String papermill = "TRUE"
        String genome = "hg38"
        
        String output_filename = "output.ipynb"
        String docker_image = ""
        Int mem_gb = 16
    }

    command {
        
        papermill $(which umap_archr.ipynb) ${output_filename} \
        -p atac_frag ${atac_frag} \
        -p genome ${genome} \
        -p papermill ${papermill}
    }

    output {
        File notebook_output = output_filename
        File archr_umap_plot = "plots/umap.png"
        File plots_zip = "plots.zip"
    }

    runtime {
        cpu : 4
        memory : mem_gb+'G'
        docker : docker_image
    }
    
    parameter_meta {
        atac_frag: {
            description: 'ATAC fragment file',
            help: 'ATAC fragments in .bedpe.gz format',
            example: 'atac.bedpe.gz'
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