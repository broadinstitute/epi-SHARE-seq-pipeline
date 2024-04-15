version 1.0

# TASK
# SHARE-qc-rna

task qc_rna {
    meta {
        version: 'v0.1'
        author: 'Mei Knudson (mknudson@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: QC RNA task'
    }

    input {
        # This function takes in input the sorted bam file produced by STARsolo
        File bam
        File mtx_tar
        Int? umi_min_cutoff = 1
        Int? gene_min_cutoff = 1
        Int? hist_min_umi = 100
        Int? hist_max_umi = 5000
        String genome_name
        String? barcode_tag = "CB"
        String? subpool
        String? prefix

        Int? cpus = 2
        Float? disk_factor = 1.0
        Float? memory_factor = 1
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_qc_rna:dev"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(bam, "G") + size(mtx_tar, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 16.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String bai = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.bam.bai"
    String barcode_metadata = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.barcode.metadata.tsv"
    String mapped_to_gene = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.reads.mapped.to.genes.txt"
    String qc_statistics = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.qc.statistics.txt"
    String umi_barcode_rank_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.umi.barcode.rank.plot.png"
    String gene_barcode_rank_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.gene.barcode.rank.plot.png"
    String gene_umi_scatter_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.gene.umi.scatter.plot.png"
    String umi_histogram_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.fragment.histogram.png"
    String monitor_log = "monitor.log"

    command <<<
        set -e

        bash $(which monitor_script.sh) | tee ~{monitor_log} 1>&2 &

        # Index bam file
        samtools index -@ ~{cpus} ~{bam} ~{bai}

        # Get STARsolo deduped statistics from mtx file
        tar -xvzf ~{mtx_tar}

        zcat matrix.mtx.gz | awk -v OFS="\t" 'NR>3{count[$2]+=$3; tot[$2]+=1}END{for (bc in count){ print bc,count[bc],tot[bc]} }' | sort -k1,1n > barcode_count_statistics_dedup.raw.tsv
        echo -e "barcode\tunique_umi\tgenes_final" > barcode_count_statistics_dedup.tsv
        awk -v pkr=~{if defined(subpool) then "_~{subpool}" else ""} -v OFS="\t" 'FNR==NR{bc[NR]=$1}FNR!=NR{print bc[$1]pkr,$2,$3; delete bc[$1]}END{for(idx in bc){print bc[idx]pkr,0,0}}' <(zcat barcodes.tsv.gz) barcode_count_statistics_dedup.raw.tsv | sort -k1,1 >> barcode_count_statistics_dedup.tsv

        # Extract barcode metadata from bam file
        python3 $(which rna_barcode_metadata.py) ~{bam} \
                                                 ~{bai} \
                                                 tmp_metadata.tsv \
                                                 ~{subpool} \
                                                 ~{"--barcode_tag " + barcode_tag}

        # Calculate FRIG (fraction of unique reads in genes)
        join -t $'\t' -e 0 -j1 <(cat tmp_metadata.tsv | (sed -u 1q;sort -k1,1)) barcode_count_statistics_dedup.tsv | \
        awk -v OFS="\t" 'NR==1{print $0,"FRIG"}NR>1{printf "%s\t%4.2f\n",$0,$9/$2}' > ~{barcode_metadata}

        # Write aggregate statistics into QC statistics file
        awk -v OFS="," 'NR>1{total+=$2; mito+=$4; unique+=$9} END {print "RNA_unique_reads_mapped_to_genes", unique; printf "RNA_FRIG,%.2f\n",unique/total; print "RNA_duplicate_reads", total-unique; printf "RNA_percent_duplicates,%.1f\n", (total-unique)/total*100; printf "RNA_percent_mitochondrial,%.1f\n", mito/total*100}' ~{barcode_metadata} > ~{qc_statistics}

        # Make QC plots
        Rscript $(which rna_qc_plots.R) ~{barcode_metadata} ~{umi_min_cutoff} ~{gene_min_cutoff} ~{hist_min_umi} ~{hist_max_umi} ~{umi_barcode_rank_plot} ~{gene_barcode_rank_plot} ~{gene_umi_scatter_plot} ~{umi_histogram_plot}
    >>>

    output {
        File rna_barcode_metadata = "~{barcode_metadata}"
        File rna_qc_statistics = "~{qc_statistics}"
        File? rna_umi_barcode_rank_plot = "~{umi_barcode_rank_plot}"
        File? rna_gene_barcode_rank_plot = "~{gene_barcode_rank_plot}"
        File? rna_gene_umi_scatter_plot = "~{gene_umi_scatter_plot}"
        File? rna_umi_histogram = "~{umi_histogram_plot}"
    }

    runtime {
        cpu : cpus
        memory : "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : "${docker_image}"
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        umi_min_cutoff: {
                description: 'UMI minimum cutoff',
                help: 'Cutoff for minimum number of UMIs required when making UMI barcode rank plot.',
                example: 10
            }
        gene_min_cutoff: {
                description: 'Gene minimum cutoff',
                help: 'Cutoff for minimum number of genes required when making gene barcode rank plot.',
                example: 10
            }
        hist_min_umi: {
                description: 'Histogram UMI minimum',
                help: 'Minimum number of UMIs per barcode when making UMI count histogram (x-axis minimum).',
                example: 100
            }
        hist_max_umi: {
                description: 'Histogram UMI maximum',
                help: 'Maximum number of UMIs per barcode when making UMI count histogram (x-axis maximum).',
                example: 5000
            }
        subpool: {
                description: 'Experiment subpool',
                help: 'Id of the sample subpool. Can be used to distinguish the 10X lanes.',
                examples: ['SS-PKR-000']
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