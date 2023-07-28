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

        pigz -d -c -p 8 ~{fragments} | awk '{size=$3-$2; if (size < 100){no_nuc+=1}; if(size >= 100 && size < 200){mono+=1}; if(size>=200){multi+=1}}END{print no_nuc*2"\n"mono*2"\n"multi*2}' > lib_sizes

        # Compute track for all insertion sizes
        bash make_track.sh -o ~{prefix}.~{genome_name}.bw ~{chrom_sizes} ~{fragments}
        # Compute track for fragments shorter than 100 nucleotides
        pigz -d -c -p 4 ~{fragments} | awk '$3-$2 < 100' | pigz -p 8 -c > no_nucleosome.bed.gz
        library_size=$(sed -n '1p;3q' $lib_sizes)
        bash make_track.sh -d $library_size -o ~{prefix}.no.nucleosome.~{genome_name}.bw ~{chrom_sizes} no_nucleosome.bed.gz
        # Compute track for mono-nuclesome fragments
        pigz -d -c -p 4 ~{fragments} | awk '$3-$2 >= 100 && $3-$2 <200' | pigz -p 8 -c > mono_nucleosome.bed.gz
        library_size=$(sed -n '2p;3q' $lib_sizes)
        bash make_track.sh -d $library_size -o ~{prefix}.mono.nucleosome.~{genome_name}.bw ~{chrom_sizes} mono_nucleosome.bed.gz
        # Compute track for fragments spanning multiple nucleosomes
        pigz -d -c -p 4 ~{fragments} | awk '$3-$2 >= 200' | pigz -p 8 -c > multi_nucleosome.bed.gz
        library_size=$(sed -n '3p;3q' $lib_sizes)
        bash make_track.sh -d $library_size -o ~{prefix}.multi.nucleosome.~{genome_name}.bw ~{chrom_sizes} multi_nucleosome.bed.gz

    >>>

    output {
        File atac_track_bigwig = "~{prefix}.~{genome_name}.bw"
        File atac_track_bigwig_no_nucleosome = "~{prefix}.no.nucleosome.~{genome_name}.bw"
        File atac_track_bigwig_mono_nucleosome = "~{prefix}.mono.nucleosome.~{genome_name}.bw"
        File atac_track_bigwig_multi_nucleosome = "~{prefix}.multi.nucleosome.~{genome_name}.bw"
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