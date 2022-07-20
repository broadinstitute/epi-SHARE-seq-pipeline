version 1.0

# TASK
# SHARE-atac-bam2bed


task share_atac_bam2bed {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC bam to bed task'
    }

    input {
        # This task takes in input the aligned bam file and rmeove the low quality reads, the extra chromosomes, marks
        # the duplicats, and convert to a bedpe file.
        Int? cpus = 8
        Int? memory_gb = 64
        File bam
        File bam_index
        File? chrom_sizes
        String genome_name
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_bam2bed"
        String? prefix


    }

    Float input_file_size_gb = size(bam, "G")
    Int mem_gb = memory_gb
    Int samtools_cpu = 6
    Int samtools_mem_gb = 8
    #Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    Int disk_gb = 200

    String filtered_chr_bam = '${default="share-seq" prefix}.filtered_chr.bam'
    String bedpe = 'tmp.bedpe'
    String final_bam = '${default="share-seq" prefix}.atac.bam2bed.alignment.cleaned.${genome_name}.bam'
    String final_bam_index = '${default="share-seq" prefix}.atac.bam2bed.alignment.cleaned.${genome_name}.bam.bai'
    String fragments = '${default="share-seq" prefix}.atac.bam2bed.fragments.${genome_name}.bgz'

    command<<<
        set -e

        # I need to do this because the bam and bai need to be in the same folder but WDL doesn't allow you to
        # co-localize them in the same path.
        ln -s ~{bam} in.bam
        ln -s ~{bam_index} in.bam.bai

        # Remove unwanted chromosomes
        chrs=$(samtools view -H in.bam | \
            grep chr | \
            cut -f2 | \
            sed 's/SN://g' | \
            grep -v chrM | \
            grep -v Y | \
            awk '{if(length($0)<6)print}')

        # Sort file by name, remove low quality reads, namesort the input bam
        samtools view -b -q 30 -f 0x2 in.bam $(echo $chrs) | \
        samtools sort -@ ~{samtools_cpu} -m ~{samtools_mem_gb}G -n -o ~{filtered_chr_bam} -

        # Convert bam to bed.gz and mark duplicates
        # Removing reads that starts and ends at the same position (duplicates)
        bedtools bamtobed -bedpe -i ~{filtered_chr_bam} | \
            sed 's/_/\t/g' | \
            awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6-5,$8}else if($10=="-"){print $1,$2-5,$6+4,$8}}' | \
            sort --parallel=~{cpus} -S 2G  -k4,4 -k1,1 -k2,2n -k3,3n | \
            uniq -c | \
            awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' > ~{bedpe}

        # Convert the bedpe file to a bam file for QC
        bedToBam -i ~{bedpe} -g ~{chrom_sizes} | \
        samtools sort -@ ~{samtools_cpu} -m ~{samtools_mem_gb}G -o ~{final_bam} -

        # and index the bam
        samtools index -@ ~{cpus} ~{final_bam}

        # Compress the bedpe file
        sort -k1,1 -k2,2n ~{bedpe} | bgzip -c > ~{fragments}

    >>>

    output {
        File atac_alignment_filtered = final_bam
        File atac_alignment_filtered_index = final_bam_index
        File atac_fragments_raw = fragments
    }

    runtime {
        cpu : cpus
        memory : mem_gb+'G'
        disks : 'local-disk ${disk_gb} SSD'
        maxRetries : 0
        docker: docker_image
    }

    parameter_meta {
        bam: {
                description: 'bam file',
                help: 'Aligned reads in bam format',
                example: 'aligned.hg38.bam'
            }
        genome_name: {
                description: 'Reference name',
                help: 'The name of the reference genome used by the aligner.',
                examples: ['hg38', 'mm10', 'both']
            }
        chrom_sizes: {
                description: 'Chromosomes size file',
                help: 'File with the length of each chromosome',
                examples: 'hg38.chrom.sizes'
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
