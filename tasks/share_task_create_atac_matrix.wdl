task create_atac_matrix {
    meta {
        version: 'v0.1'
        author: 'Siddarth Wekhande (swekhand@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: create ATAC count matrix task'
    }
    
    input {
        # This takes in ATAC fragment file and a peak set, and outputs a h5 count matrix
        File atac_fragment_file
        File peak_file
        File output_h5
        
        String docker_image  
        Int mem_gb = 8
    }
    
    
    
    command {
        set -e Rscript $(which create_atac_matrix.R) ${atac_fragment_file} ${peak_file} ${output_h5}
    }
    
    output {
        File count_matrix_h5= glob('*.h5')
    }

    runtime {
        memory : '${mem_gb} GB'
        docker: docker_image
    }
    
    parameter_meta {
        atac_fragment_file: {
                description: 'ATAC fragment file',
                help: 'ATAC fragments stores as a bed/bedpe file',
                example: 'atac.bedpe'
            },
        peak_file: {
                description: 'Consensus peak set',
                help: 'Consensus peak set, bed file (mostly cCRE bed file)'
                example: 'ccre.bed'
            },
        output_h5: {
                description: 'Output h5',
                help: 'Name of the output h5 file'
                example: 'output.h5'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: R and R libraries [SummarizedExperiment, dplyr, GenomicRanges, Matrix, DropletUtils]'
                example: ['put link to gcr or dockerhub']
            }
    }
}