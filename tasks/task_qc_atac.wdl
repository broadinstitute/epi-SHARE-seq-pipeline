version 1.0

# TASK
# qc-atac

task qc_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard IGVF pipeline: ATAC qc statistics task'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File fragment_file
        File? fragment_file_index
        File chrom_sizes
        File gtf
        # File? tss_bed
        Int fragment_min_cutoff = 10
        String? prefix

        # Runtime
        Int cpus = 10
        Float disk_factor = 10.0
        Float memory_factor = 0.3
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_atac:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragment_file, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 32.0 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(100.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    # pdf string needed as required input to Picard CollectInsertSizeMetrics
    String fragment_size_distribution_plot = "~{prefix}_fragment_size_distribution.png"
    String tss_enrichment_library_plot = "~{prefix}_TSS_enrichment.png"
    String fraction_of_duplicates_distribution_plot = "~{prefix}_fraction_of_duplicates_distribution.png"
    String fraction_of_mito_distribution_plot = "~{prefix}_fraction_of_mitochondrial_fragments_distribution.png"
    String library_tss_overlap = "~{prefix}_frac_overlap_TSS.txt"
    String library_tsse = "~{prefix}_library_TSS.txt"
    String knee_plot = "~{prefix}_knee_plot.png"
    String n_fragment_vs_tss_enrichment_plot = "~{prefix}_n_fragment_vs_TSS_enrichment.png"
    String n_fragment_vs_tss_enrichment_filtered_plot = "~{prefix}_n_fragment_vs_TSS_enrichment_filtered.png"
    String umap_leiden_plot = "~{prefix}_umap_leiden.png"
    String snapatac2_h5ad = "~{prefix}_snap.h5ad"


    command <<<
        qc_atac \
            --fragment_file ~{fragment_file} \
            --chrom_sizes ~{chrom_sizes} \
            --compressed_gtf_file ~{gtf} \
            --min_frag_cutoff ~{fragment_min_cutoff} \
            --prefix ~{prefix}
    >>>

    output {
        File atac_fragment_size_distribution_plot = fragment_size_distribution_plot
        File atac_tss_enrichment_library_plot = tss_enrichment_library_plot
        File atac_fraction_of_duplicates_distribution_plot = fraction_of_duplicates_distribution_plot
        File atac_fraction_of_mito_distribution_plot = fraction_of_mito_distribution_plot
        File atac_knee_plot = knee_plot
        File atac_n_fragment_vs_tss_enrichment_plot = n_fragment_vs_tss_enrichment_plot
        File atac_n_fragment_vs_tss_enrichment_filtered_plot = n_fragment_vs_tss_enrichment_filtered_plot
        File atac_umap_leiden_plot = umap_leiden_plot
        File atac_snapatac2_h5ad = snapatac2_h5ad
        File atac_barcode_metrics = "~{prefix}_barcode_metrics.csv"
        Float atac_library_tss_overlap = read_float(library_tss_overlap)
        Float atac_library_tsse = read_float(library_tsse)
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2',
                examples: '4'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
        fragment_file: {
            description: "Fragment file",
            help: "The input fragment file containing the raw sequencing reads.",
            example: "sample.fragments.tsv.gz"
        }
        fragment_file_index: {
            description: "Fragment file index",
            help: "Optional index file for the fragment file.",
            example: "sample.fragments.tsv.gz.tbi"
        }
        chrom_sizes: {
            description: "Chromosome sizes file",
            help: "File containing the sizes of the chromosomes.",
            example: "hg38.chrom.sizes"
        }
        gtf: {
            description: "GTF file",
            help: "GTF file containing gene annotations.",
            example: "genes.gtf"
        }
        prefix: {
            description: "Output prefix",
            help: "Prefix for the output files.",
            example: "sample"
        }
    }
}
