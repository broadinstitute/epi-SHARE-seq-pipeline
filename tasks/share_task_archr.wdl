version 1.0

task archr {
    meta {
        version: 'v0.1'
        author: 'Kevin Dong (kdong@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: run archr task'
    }

    input {
        #This tasks takes in an ATAC fragment file, processes using ArchR and creates plots

        #ArchR parameters
        File atac_frag
        String genome

        #ArchR QC
        Int min_tss = 4
        Int min_frags = 100
        String add_tile_mat = "TRUE"
        String add_gene_score_mat = "TRUE"

        #ArchR Doublet paramaters
        Int doublet_k = 10
        String doublet_knn_method = "UMAP"
        Int lsi_method = 1

        String copy_arrow_files = "TRUE"
        String iter_LSI_matrix = "TileMatrix"
        Int threads = 1
        String prefix = "prefix"

        #ArchR Plots parameters
        String marker_features_test = "wilcoxon"
        String heatmap_transpose = "TRUE"
        Int heatmap_label_n = 5

        String heatmap_cutoff = "'FDR <= 0.01 & Log2FC >= 0.5'" # Fix - use two float as parameters and create the string inside the task

        #papermill specific parameters
        String papermill = "TRUE"

        String output_filename = "output.ipynb"
        String docker_image = "swekhande/dorcs:epi-umap-archr"
        Int? mem_gb = 32
        Int? disk_gb = 50
    }

    #Plot filepaths
    String doublet_summary_plot = 'QualityControl/${prefix}/${prefix}-Doublet-Summary.pdf'
    String fragment_size_dist_plot = 'QualityControl/${prefix}/${prefix}-Fragment_Size_Distribution.pdf'
    String TSS_uniq_frags_plot = 'QualityControl/${prefix}/${prefix}-TSS_by_Unique_Frags.pdf'
    String TSS_uniq_frags_filtered_plot = 'plots/${prefix}.atac.archr.TSS_fragment_qc.${genome}.png'
    String heatmap_plot = 'plots/${prefix}.atac.archr.heatmap.${genome}.png'
    String umap_plot = 'plots/${prefix}.atac.archr.umap.${genome}.png'

    #Other filepaths
    String arrow_file = '${prefix}.arrow'
    String archr_rds = '${prefix}.atac.archr.rds.${genome}.rds'
    String plots_zip_dir = 'plots.zip'


    command {

        papermill $(which archr_notebook.ipynb) ${output_filename} \
        -p atac_frag ${atac_frag} \
        -p genome ${genome} \
        -p papermill ${papermill} \
        -p min_tss ${min_tss} \
        -p min_frags ${min_frags} \
        -p add_tile_mat ${add_tile_mat} \
        -p add_gene_score_mat ${add_gene_score_mat} \
        -p doublet_k ${doublet_k} \
        -p doublet_knn_method ${doublet_knn_method} \
        -p lsi_method ${lsi_method} \
        -p copy_arrow_files ${copy_arrow_files} \
        -p iter_LSI_matrix ${iter_LSI_matrix} \
        -p threads ${threads} \
        -p prefix ${prefix} \
        -p marker_features_test ${marker_features_test} \
        -p heatmap_transpose ${heatmap_transpose} \
        -p heatmap_label_n ${heatmap_label_n} \
        -p heatmap_cutoff ${heatmap_cutoff}
    }

    output {
        File notebook_output = output_filename
        #File archr_umap_plot = umap_plot
        #File? archr_heatmap_plot = heatmap_plot
        #File archr_TSS_uniq_frags_plot = TSS_uniq_frags_plot
        #File archr_TSS_uniq_frags_filtered_plot = TSS_uniq_frags_filtered_plot
        #File archr_fragment_size_dist_plot = fragment_size_dist_plot
        #File archr_doublet_plot = doublet_summary_plot

        #File plots_zip = plots_zip_dir
        #File archr_arrow = arrow_file
        #File archr_obj = archr_rds
    }

    runtime {
        cpu : 4
        memory : mem_gb+'G'
        docker : docker_image
        disks : 'local-disk ${disk_gb} LOCAL'
        maxRetries : 0
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
            help: 'Docker image for preprocessing step.' ,
            example: ['put link to gcr or dockerhub']
        }

        min_tss: {
            description: 'Min TSS enrichment score',
            help: 'The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering',
            example: 4
        }

        min_frags: {
            description: 'Min number of mapped fragments',
            help: 'The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use',
            example: 1000
        }

        add_tile_mat: {
            description: 'Compute Tile Matrix if TRUE',
            help: 'A boolean value indicating whether to add a "Tile Matrix" to each ArrowFile.',
            example: 'TRUE'
        }

        add_gene_score_mat: {
            description: 'Compute Gene Score Matrix if TRUE',
            help: 'A boolean value indicating whether to add a Gene-Score Matrix to each ArrowFile.',
            example: 'TRUE'
        }

        doublet_k: {
            description: 'Number of simulated neighbors to be considered doublet',
            help: 'The number of cells neighboring a simulated doublet to be considered as putative doublets.',
            example: 10
        }

        doublet_knn_method: {
            description: 'Embedding method for doublet detection',
            help: 'Refers to the embedding to use for nearest neighbor search.',
            examples: ['UMAP','LSI']
        }

        lsi_method: {
            description: 'Order of operations of TF-IDF normalization (see ArchR manual)',
            help: 'A number or string indicating the order of operations in the TF-IDF normalization. Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf"',
            examples: [1,2,3]
        }

        copy_arrow_files: {
            description: 'Makes a copy of arrow files',
            help: 'Save a copy of arrow files in the ArchR project (recommended)',
            example: 'TRUE'
        }

        iter_LSI_matrix: {
            description: 'Data matrix to retrieve',
            help: 'The name of the data matrix to retrieve from the ArrowFiles associated with the ArchRProject.',
            examples: ['PeakMatrix','TileMatrix']
        }

        threads: {
            description: 'Number of threads to run ArchR',
            help: 'Set threads to run ArchR. For now, recommended to run on single (1) thread.',
            example: 1
        }

        prefix: {
            description: 'Project name',
            help: 'String used to name your project and associated file names',
            example: "shareseq"
        }

        marker_features_test: {
            description: 'Pairwise test method',
            help: 'The name of the pairwise test method to use in comparing cell groupings to the null cell grouping during marker feature identification.',
            examples: ['wilcoxon', 'ttest', 'binomial']
        }

        heatmap_transpose: {
            description: 'Boolean to transpose heatmap',
            help: 'Plots genes on columns if TRUE',
            example: 'TRUE'
        }

        heatmap_label_n: {
            description: 'Top n genes to label',
            help: 'Extracts the top n upregulated genes in cluster to label on heatmap',
            example: 5
        }

        heatmap_cutoff: {
            description: 'Cut-off applied to genes in heatmap',
            help: 'Cut-off has to be specified in string format',
            example: 'FDR <= 0.01 & Log2FC >= 0.5'
        }
    }
}
