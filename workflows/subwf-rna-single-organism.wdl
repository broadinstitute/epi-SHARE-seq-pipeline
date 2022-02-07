version 1.0


# Import the tasks called by the pipeline
import "../tasks/share_task_star.wdl" as share_task_align
import "../tasks/share_task_update_rgid.wdl" as share_task_update_rgid
import "../tasks/share_task_count_rna.wdl" as share_task_feature_counts
import "../tasks/share_task_group_umi.wdl" as share_task_group_umi

workflow wf_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) and Sai Ma @ Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: Sub-workflow to process the RNA portion of SHARE-seq libraries.'
    }

    input {
        # RNA Sub-worflow inputs

        # Align
        File read1
        File idx_tar
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16
        String docker = "polumechanos/share-seq"
        # Update RGID
        Boolean multimappers = false
        # Assign features
        Boolean include_multimappers = false
        Boolean include_introns = false
        File gtf
        String gene_naming = "gene_name"
        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff = 100
        # Lib_size QC
        Boolean qc = false
        File genes_annotation_bed
    }

    call share_task_align.share_rna_align as align {
        input:
            fastq_R1 = read1,
            genome_name = genome_name,
            genome_index_tar = idx_tar,
            prefix = prefix,
            cpus = cpus
    }

    call share_task_update_rgid.share_rna_update_rgid as update_rgid{
        input:
            bam = align.rna_alignment,
            multimapper = multimappers,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_feature_counts.feature_counts_rna as count{
        input:
            multimapper = include_multimappers,
            intron = include_introns,
            bam = update_rgid.rna_reheaded_alignment,
            gtf = gtf,
            gene_naming = gene_naming,
            genome_name = genome_name,
            prefix = prefix
    }

    call share_task_group_umi.group_umi_rna as group_umi{
        input:
            bam = count.rna_featurecount_alignment,
            mode = mode,
            cutoff = cutoff,
            remove_single_umi = remove_single_umi,
            genome_name = genome_name,
            prefix = prefix
    }

    # TODO: the genes annotation needs to be passed just like the gtf
    # Check which one is the bam in input
    call qc_libsize_rna{
        input:
            qc = qc,
            bam = count.rna_featurecount_alignment,
            genes_annotations_bed = genes_annotation_bed,
            genome_name = genome_name,
            prefix = prefix
    }

    call qc_lib_size {
        input:
            raw_counts = group_umi.rna_umi_counts_unfiltered,
            filtered_counts = group_umi.rna_umi_counts_filtered,
            cutoff = cutoff,
            genome_name = genome_name,
            prefix = prefix
    }

    call rna_gene_cell_matrix {
        input:
            filtered_bed = group_umi.rna_umi_bed_filtered,
    }

    output {
        File share_rna_alignment_raw = align.rna_alignment
        File share_rna_alignment_index = align.rna_alignment_index
        File share_rna_alignment_log = align.rna_alignment_log

        File share_rna_reheaded_alignment = update_rgid.rna_reheaded_alignment
        File share_rna_reheaded_alignment_index = update_rgid.rna_reheaded_alignment_index

        File share_rna_featurecount_alignment = count.rna_featurecount_alignment
        File share_rna_featurecount_alignment_index = count.rna_featurecount_alignment_index
        File share_rna_featurecount_log = count.rna_featurecount_log
        File share_rna_featurecount_txt = count.rna_featurecount_txt
        File share_rna_featurecount_summary = count.rna_featurecount_summary

        File share_rna_umi_barcodes = group_umi.rna_umi_barcodes_filtered
        File share_rna_umi_bed_filtered = group_umi.rna_umi_bed_filtered
        File share_rna_umi_bed_unfiltered = group_umi.rna_umi_bed_unfiltered
        File share_rna_umi_counts_filtered = group_umi.rna_umi_counts_filtered
        File share_rna_umi_counts_unfiltered = group_umi.rna_umi_counts_unfiltered

        File rna_read_distribution_txt = qc_libsize_rna.read_distribution_rna
        File rna_read_distribution2_txt = qc_libsize_rna.read_distribution2_rna
        File rna_read_distribution_plot = qc_libsize_rna.read_distribution_plot_rna

        File rna_lib_size_count = qc_lib_size.lib_size_counts
        File rna_duplicates_log = qc_lib_size.lib_size_log
        Array[File] rna_lib_size_plots = qc_lib_size.plots

        File h5_matrix = rna_gene_cell_matrix.h5_matrix
        Array[File] umi_gene_qc_plots = rna_gene_cell_matrix.umi_qc_plots
    }
}

# TASKS










task qc_libsize_rna {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: group umi rna task'
    }

    input {
        # This function takes in input the bam file produced by the STAR
        # aligner and group UMI reads whilst filtering duplicates.
        Boolean qc = false
        File bam
        File genes_annotations_bed
        String genome_name
        String? prefix
        String docker_image
    }

    #Float input_file_size_gb = size(input[0], "G")
    Int mem_gb = 8
    Int disk_gb = 50
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)

    command {
        set -e
        # Calculate gene body coverage and reads distribution
        mv ${bam} in.bam
        INPUT=in.bam

        mv ${genes_annotations_bed} gene.bed.gz

        gunzip gene.bed.gz

        samtools index in.bam

        if [[ '${qc}' == 'true' ]]; then
            $(which samtools) view -s 0.01 -o temp.1pct.bam ${bam}
            INPUT="temp.1pct.bam"
        fi

        ## TODO: Where is this python script?
        # Calculate read distribution
        python3 $(which read_distribution.py) -i $INPUT -r gene.bed > ${prefix + "."}rna.${genome_name}.reads_distribution.txt 2>>./Run.log

        # The two files created here are necessary for the
        # Read_distribution.R script.
        tail -n +5 ${prefix + "."}rna.${genome_name}.reads_distribution.txt | head -n -1 > temp1.txt
        head -n 3  ${prefix + "."}rna.${genome_name}.reads_distribution.txt | grep Assigned | sed 's/Total Assigned Tags//g' | sed 's/ //g' > temp2.txt

        # Plot reads distribution
        Rscript $(which Read_distribution.R) . ${prefix + "."}rna.${genome_name} --save

    }

    output {
        File read_distribution_rna= "${prefix + "."}rna.${genome_name}.reads_distribution.txt"
        File read_distribution2_rna= "${prefix + "."}rna.${genome_name}.reads_distribution2.txt"
        File read_distribution_plot_rna= "${prefix + "."}rna.${genome_name}.reads_distribution.pdf"
    }

    runtime {
        #cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        qc: {
                description: 'QC run flag',
                help: 'Flag to set if you are doing a qc run.',
                default: false,
                example: [true, false]
            }
        genes_annotations_bed: {
                description: 'Genes annotations in bed format',
                help: 'The genes annotations to use when calculating the reads distributions.',
                example: 'hg38.UCSC_RefSeq.bed'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name genome reference used to align.',
                example: ['hg38', 'mm10', 'hg19', 'mm9']
            }
        prefix: {
                description: 'Prefix for output files',
                help: 'Prefix that will be used to name the output files',
                example: 'MyExperiment'
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: samtools',
                example: ['put link to gcr or dockerhub']
            }
    }
}


task qc_lib_size {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: RNA library size task'
    }

    input {
        # This task computs the the library size for the library.

        File raw_counts
        File filtered_counts
        Int cutoff
        String genome_name
        String? prefix
        String docker_image


    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    Float input_file_size_gb = size(filtered_counts, "G")
    Int mem_gb = 16


    command {
        set -e

        # TODO remove the hard coded file paths from R scripts
        # TODO create only one R script that uses the parameters to discriminate
        # Estimate lib size
        # both
        #Rscript $(which lib_size_sc_V5_species_mixing.R)./ '${prefix + '.'}atac.${genome_name}' ${cutoff} atac --save
        # hg38/mm10
        Rscript $(which lib_size_sc_V5_single_species.R) ${raw_counts} ${filtered_counts} ${cutoff} ${genome_name} ATAC --save

    }

    output {
        File lib_size_counts = glob('*.libsize.counts.csv')[0]
        File lib_size_log = glob('*.dups.log')[0]
        Array[File] plots = glob('*.pdf')
    }

    runtime {
        #cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        raw_counts: {
                description: 'Barcode count csv',
                help: 'Barcode counts from raw bam in csv format.',
                example: 'raw.counts.csv'
            }
        filtered_counts: {
            description: 'Barcode count csv',
            help: 'Barcode counts from filtered bam in csv format.',
            example: 'filtered.counts.csv'
        }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }


}

task rna_gene_cell_matrix {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: RNA gene cell matrix'
    }

    input {
        # This task computs the the gene by barcode matrix.

        File filtered_bed
        String docker_image


    }

    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 50
    Float input_file_size_gb = size(filtered_bed, "G")
    Int mem_gb = 16


    command {
        set -e
        Rscript $(which UMI_gene_perCell_plot_v2.R) ${filtered_bed} --save

    }

    output {
        File h5_matrix = glob('*.h5')[0]
        Array[File] umi_qc_plots = glob('*.pdf')
    }

    runtime {
        #cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
            filtered_bed: {
            description: 'Barcode count csv',
            help: 'Barcode counts from filtered bam in csv format.',
            example: 'filtered.counts.csv'
        }
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz',
                example: ['put link to gcr or dockerhub']
            }
    }
}
