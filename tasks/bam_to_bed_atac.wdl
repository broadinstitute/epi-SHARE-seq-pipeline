task bam_to_bed_atac {
    meta {
        version: 'v0.1'
        author: 'Eugenio Mattei (emattei@broadinstitute.org) at Broad Institute of MIT and Harvard'
        description: 'Broad Institute of MIT and Harvard SHARE-Seq pipeline: ATAC bam to bed task'
    }
    
    input {
        # This function takes in input the raw fastqs from Novaseq or
        # Nextseq and perform trimming, adapter removal, appending
        # cell barcode to the read name, and splitting ATAC and RNA.
        
        File bam
        String genome_name
        String type= "atac"
        File chrom_sizes
        String docker_image
        Int cpus= 4
        
    }
    
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 8.0
    Int disk_gb = round(20.0 + 4 * input_file_size_gb)
    
    String filtered_chr_bam= 'filtered_chr.bam'
    String bedpe= 'tmp.bedpe'
    String final_bam= '${type}.cleaned.${genome_name}.bam'
    String final_bam_index= '${type}.cleaned.${genome_name}.bam.bai'
    String final_bedpe= '${type}.cleaned.${genome_name}.bedpe.gz'

    command {
        set -e
        
        # Remove unwanted chromosomes
        chrs=`samtools view -H ${bam} | \
            grep chr | \
            cut -f2 | \
            sed 's/SN://g' | \
            grep -v chrM | \
            grep -v Y | \
            awk '{if(length($0)<6)print}'`
            
        # Sort file by name, remove low quality reads, namesort the input bam
        samtools view \
            -b \
            -q 30 \
            -f 0x2 ${bam}  `echo $chrs` | \
        samtools sort -@ ${cpus} -m ${mem_gb}G -n -o ${filtered_chr_bam}

        # Convert bam to bed.gz and mark duplicates
        bedtools bamtobed -i ${filtered_chr_bam} \
            -bedpe | \
        sed 's/_/\t/g' | \
        awk -v OFS="\t" '{if($10=="+"){print $1,$2+4,$6+4,$8}else if($10=="-"){print $1,$2-5,$6-5,$8}}' |\
        sort --parallel=${cpus} -S ${mem_gb}G  -k4,4 -k1,1 -k2,2 -k3,3 | \
        uniq -c | \
        awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' > ${bedpe}
        
        
        # Convert to a bam file for QC
        bedToBam -i ${bedpe} -g ${chrom_sizes} | \
        samtools sort -@ ${cpus} -m ${mem_gb}G - > ${final_bam}
        
        # and index the bam
        samtools index -@ ${cpus} ${final_bam}
        
        # Compress the bedpe file
        pigz --fast -c -p ${cpus} ${bedpe} > ${final_bedpe}
        
    }
    
    output {
        File bam_filtered= final_bam
        File bam_filtered_index= final_bam_index
        File bedpe_cleaned= final_bedpe
    }

    runtime {
        cpu : ${cpus}
        memory : '${mem_gb} GB'
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0 
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
                help: 'The name of the reference genome used by the aligner.'
                examples: ['hg38', 'mm10', 'both']
            },
        chrom_sizes: {
                description: 'Chromosomes size file',
                help: 'File with the length of each chromosome'
                examples: 'hg38.chrom.sizes'
            },
        cpus: {
                description: 'Number of cpus',
                help: 'Set the number of cpus useb by bowtie2'
                examples: '4'
            },
        docker_image: {
                description: 'Docker image.',
                help: 'Docker image for preprocessing step. Dependencies: python3 -m pip install Levenshtein pyyaml Bio; apt install pigz'
                example: ['put link to gcr or dockerhub']
            }
    }


}
