version 1.0

# TASK
# Make bigwig track from fragment file

task make_track {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC make track task'
    }

    input {
        # This function takes in input the raw and filtered bams
        # and compute some alignment metrics along with the TSS
        # enrichment plot.
        File? chrom_sizes
        File? fragments
        String? genome_name
        String? prefix

        # Runtime
        Int? cpus = 8
        Float? disk_factor = 4
        Float? memory_factor = 0.3
        String docker_image = "us.gcr.io/buenrostro-share-seq/task_make_track"
    }

    # Determine the size of the input
    Float input_file_size_gb = size(fragments, "G")

    # Determining memory size base on the size of the input files.
    Float mem_gb = 8 + memory_factor * input_file_size_gb

    # Determining disk size base on the size of the input files.
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"


    String monitor_log = "atac_qc_monitor.log"


    command<<<
        set -e

        # I am not writing to a file anymore because Google keeps track of it automatically.
        bash $(which monitor_script.sh) 1>&2 &

        pigz -c -d -p 8 ~{fragments} | \
        awk -v OFS="\t" '{print $1,$2-4,$2+4"\n"$1,$3-4,$3+4}' | \
        sort --parallel=4 -k1,1 -k2,2n > tn5_insertions.bed

        insertion_number=$(wc -l tn5_insertions.bed)
        scale_factor=$(bc <<< "scale=6;10000000/$(echo $insertion_number)")

        bedtools merge -i insertions.bed -c 1 -o count | \
        awk -v scaling=$scale_factor -v OFS="\t" '{$4=$4*scaling; print $0}' | \
        sort --parallel=8 -k1,1 -k2,2n > ~{prefix}.bedGraph

        bedGraphToBigWig ~{prefix}.bedGraph ~{chrom_sizes} ~{prefix}.~{genome_name}.bw

    >>>

    output {
        File atac_track_bigwig = "~{prefix}.bigwig"
    }

    runtime {
        cpu: cpus
        disks: "local-disk ${disk_gb} ${disk_type}"
        docker: "${docker_image}"
        maxRetries: 0
        memory: "${mem_gb} GB"
    }

    parameter_meta {
        fragments: {
                description: 'Fragment file',
                help: 'Fragment file.',
                example: 'fragments.tsv.gz'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
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
    }
}
