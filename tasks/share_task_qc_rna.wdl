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
        Int? cutoff
        String genome_name
        String? prefix
        
        Int? cpus = 16
        Float? disk_factor = 8.0
        Float? memory_factor = 0.15

        #String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_qc_rna"
        String docker_image = "mknudson/share_task_qc_rna:test"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(bam, "G")

    # Determining memory size based on the size of the input files.
    Float mem_gb = 5.0 + memory_factor * input_file_size_gb

    # Determining disk size based on the size of the input files.
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # Determining disk type based on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    String assay = "RNA"
    String bai = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.bam.bai"
    String barcode_metadata = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.barcode.metadata.txt"
    String duplicates_log = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.duplicates.log.txt"
    String barcode_rank_plot = "~{default="share-seq" prefix}.qc.rna.~{genome_name}.barcode.rank.plot.png"

    command <<<
        set -e

        bash $(which monitor_script.sh) > monitoring.log &

        # Index bam file
        samtools index -@ ~{cpus} ~{bam} ~{bai} 
         
        # Extract barcode metadata (total counts, unique counts, duplicate counts, genes, percent mitochondrial) from bam file
        python3 $(which rna_barcode_metadata.py) ~{bam} ~{bai} ~{default="share-seq" prefix} ~{barcode_metadata}
    
        # Make duplicates log from barcode metadata file
        awk '{total+=$2; duplicate+=$3; unique+=$4} END {print "total reads:", total; print "unique reads:", unique; print "duplicate reads:", duplicate}' ~{barcode_metadata} > ~{duplicates_log}

        # Get UMI column of barcode_metadata file, remove header, sort
        cut -f4 ~{barcode_metadata} | tail -n +2 | sort -rn > umis_per_barcode
        # Make barcode rank plot 
        Rscript $(which barcode_rank_plot.R) umis_per_barcode ~{cutoff} ~{genome_name} ~{assay} ~{barcode_rank_plot}
    >>>

    output {
        File monitor_log = "monitoring.log"
        File rna_barcode_metadata = "${barcode_metadata}"
        File rna_duplicates_log = "${duplicates_log}"
        File rna_barcode_rank_plot = "${barcode_rank_plot}"   
    }

    runtime {
        cpu : cpus
        memory : "~{mem_gb} GB"
        memory_retry_multiplier: 2
        disks: "local-disk ~{disk_gb} ~{disk_type}"
        docker : docker_image
        maxRetries:1
    }

    parameter_meta {
        bam: {
                description: 'Alignment bam file',
                help: 'Aligned reads in bam format.',
                example: 'hg38.aligned.bam'
            }
        cutoff: {
                description: 'UMI cutoff',
                help: 'Cutoff for number of UMIs required when plotting barcode statistics.',
                example: 100
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
