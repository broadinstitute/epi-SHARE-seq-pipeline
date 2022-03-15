version 1.0

# Import the sub-workflow for preprocessing the fastqs.
import "workflows/subwf-atac-single-organism.wdl" as share_atac
import "workflows/subwf-rna-single-organism.wdl" as share_rna
import "workflows/subwf-find-dorcs.wdl" as find_dorcs

# WDL workflow for SHARE-seq

workflow ShareSeq {

    input {
        # Common inputs
        String prefix = "shareseq-project"
        String genome_name
        Int? cpus = 16


        # ATAC specific inputs
        File chrom_sizes
        File read1_atac
        File read2_atac
        File idx_tar_atac
        File tss_bed
        Int? cpus_atac
        Int cutoff_atac

        # RNA specific inputs
        Boolean multimappers = false
        Boolean include_multimappers = false
        Boolean include_introns = false
        File read1_rna
        File genes_annotation_bed
        File gtf
        File idx_tar_rna
        Int? cpus_rna
        String? gene_naming = "gene_name"

        # Group UMI
        Boolean remove_single_umi = true
        String mode = "fast"
        Int cutoff_rna = 100

        # Lib_size QC
        Boolean qc = false


        # DORCs specific inputs
        File peak_set
        Int? cpus_dorcs = 4
        String save_plots_to_dir = "TRUE"
        String output_filename = "output.ipynb"

        # Seurat filters
        Int minFeature_RNA = 200
        Int maxFeature_RNA = 2500
        Float percentMT_RNA = 5
        Int minCells_RNA = 3

        # DORCs filter
        Int dorcGeneCutOff = 10
        Float fripCutOff = 0.3
        Float corrPVal = 0.05
        Int topNGene = 20

        # Regulatory region around TSS. Default is +/- 50Kb
        Int windowPadSize = 50000
        Int bootstraps = 100

        String docker_image_dorcs = "polumechanos/dorcs_task_find_dorcs"
        Int mem_gb_dorcs = 16
    }

    call share_rna.wf_rna as rna{
        input:
            read1 = read1_rna,
            idx_tar = idx_tar_rna,
            prefix = prefix,
            genome_name = genome_name,
            cpus = cpus_rna,
            # Update RGID
            multimappers = multimappers,
            # Assign features
            include_multimappers = include_multimappers,
            include_introns = include_introns,
            gtf = gtf,
            gene_naming = gene_naming,
            # Group UMI
            remove_single_umi = remove_single_umi,
            mode = mode,
            cutoff = cutoff_rna,
            # Lib_size QC
            qc = qc,
            genes_annotation_bed = genes_annotation_bed
    }
    call share_atac.wf_atac as atac{
        input:
            read1 = read1_atac,
            read2 = read2_atac,
            chrom_sizes = chrom_sizes,
            idx_tar = idx_tar_atac,
            tss_bed = tss_bed,
            prefix = prefix,
            genome_name = genome_name,
            cutoff = cutoff_atac,
            cpus = cpus_atac
    }
    call find_dorcs.wf_dorcs as dorcs{
        input:
            rna_matrix = rna.share_rna_h5_matrix,
            atac_fragments = atac.share_atac_fragments_filtered,
            peak_file = peak_set,

            genome = genome_name,
            n_cores = cpus_dorcs,
            save_plots_to_dir = save_plots_to_dir,
            output_filename = output_filename,

            minFeature_RNA = minFeature_RNA,
            maxFeature_RNA = maxFeature_RNA,
            percentMT_RNA = percentMT_RNA,
            minCells_RNA = minCells_RNA,

            dorcGeneCutOff = dorcGeneCutOff,
            fripCutOff = fripCutOff,
            corrPVal = corrPVal,
            topNGene = topNGene,

            # Regulatory region around TSS. Default is +/- 50Kb
            windowPadSize = windowPadSize,
            bootstraps = bootstraps,
            mem_gb = mem_gb_dorcs
    }

}


task merge_fastqs{
    input{
        Array[File] atac_read1
        Array[File] atac_read2
        Array[File] rna_read1
        String? prefix
    }
    command{
        cat ${sep=' ' atac_read1} > ${prefix + "."}merged.atac.R1.fq.gz
        cat ${sep=' ' atac_read2} > ${prefix + "."}merged.atac.R2.fq.gz
        cat ${sep=' ' rna_read1} > ${prefix + "."}merged.rna.R1.fq.gz
    }
    output{
        File merged_atac_fastq_R1 = glob('*.merged.atac.R1.fq.gz')[0]
        File merged_atac_fastq_R2 = glob('*.merged.atac.R2.fq.gz')[0]
        File merged_rna_fastq_R1 = glob('*.merged.rna.R1.fq.gz')[0]
    }
}


# Task to report errors to user.
# From https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl
task raise_exception {
  input {
    String msg
    Array[String]? vals
  }
  command {
    echo -e "\n* Error: ${msg}\n" >&2
    echo -e "* Vals: ${sep=',' vals}\n" >&2
    exit 2
  }
  output {
    String error_msg = '${msg}'
  }
  runtime {
    maxRetries : 0
    cpu : 1
    memory : '2 GB'
    time : 1
    disks : 'local-disk 10 SSD'
  }
}
