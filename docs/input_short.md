# Input JSON

The input JSON file must include all the information needed to run this pipeline. Hence, it must include the absolute paths to all the control and experimental FASTQ files, paths to all the genomic data files needed for this pipeline, as well as the parameters and metadata needed for running this pipeline. If the parameters are not specified in the input JSON file, default values will be used. We provide a set of template JSON files: [minimal](../example_input_json/input-short-share.json) and [full](../example_input_json/input-complete-share.json). We recommend using a minimal template rather than a full template. A full template includes all parameters of the pipeline with default values defined.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

# IGVF jamboree
The input JSONs used to run the SHARE-seq and 10x samples can be found at the following locations:
1) SHARE-seq human [syn...](putlinktosynapse)
2) SHARE-seq mouse [syn...](putlinktosynapse)
3) 10X Multiome human [syn...](putlinktosynapse)
3) 10X Multiome mouse [syn....](putlinktosynapse)

# Checklist

Mandatory parameters.

1) Modality
    * `share.pipeline_modality`:
        * `full`: the **default**. Will run the entire pipeline.
        * `count_only`: Will stop after generating the `fragment` file and/or the `barcode x gene` matrix.
        * `no_align`: Will perform pre-processing of the raw FASTQs and stop. 
    > **IMPORTANT**: Because a bug in ArchR the `full` mode will not work with `mm39`

2) Chemistry
    * `share.chemistry`:`shareseq` for SHARE-seq, `10x_multiome` for 10x multiome and `10x_v2` for not multiomic 10x data.


3) Reference genome tsv
    * `share.genome_tsv`: Define the TSV file containing the genome annotations necessary to run the pipeline. 
    * The **bolded** options in the table below represent the defaults for the two organisms and the different GENCODE and assembly versions.

        Genome|URL
        -|-
        **hg38_gencode_v43**|`gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human_v43/Homo_sapiens_genome_files_hg38_v43.tsv`
        hg38_gencode_v42|`gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human/GRCh38_genome_files_hg38.tsv`
        **mm10**|`gs://broad-buenrostro-pipeline-genome-annotations/mm10/mm10_genome_files_STARsolo.tsv`
        mm39_IGVF_gencode_v32|`gs://broad-buenrostro-pipeline-genome-annotations/IGVF_mouse_v32/Mus_musculus_genome_files_mm39_v32.tsv`

        > **IMPORTANT**: Because a bug in ArchR the `full` mode will not work with `mm39` and `mm39_IGVF_gencode_v32`

    * To build a new TSV file from use your own FASTA (`.fa` and `.2bit`), see [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/build_genome_database.md).
    * You can also define genome-specific parameters (defined in a genome TSV file) in an input JSON file. Parameters defined in an input JSON file will override those defined in a genome TSV file. For example, you can simply replace a blacklist file while keeping all other parameters in a genome TSV file for `hg38`.
        ```javascript
        {
            "share.tss_bed": "/new/genome/data/new_tss.bed.gz",
            "share.peak_set": "File? (optional)"
        }
        ```

4) [Input files](#input-files)
    * See [this](#input-files) for how to define FASTQs for your sample.

5) Important parameters
    * If you are running samples generated from `cellular/nuclear subpools` it is recommended to define the following parameters.
        * `share.pkr`: Cellular/nuclear subpool identifier. The pipeline will add this parameter to the cell barcode in order to prevent mixing barcodes from different subpools.
        >IMPORTANT: Processing `subpools` at the moment is a two-step process. Each `subpool` needs to be processed independently with the `share.pkr` to track the `subpools` provenance and is only aggregated afterwards.
    * `share.prefix`: Prefix used to name the output files.
    * `share.align_multimappers`: This is an integer indicating how many multimappers are allowed when aligning gDNA using bowtie2. The default is to not use multimappers,k but it is recommended to use a value of `5` as recommended by the ENCODE DAC.
    * `share.barcode_tag_fragments`: When processing samples from different `subpools`, it is important to set this value to `XC`. This is the defaut for the `shareseq` chemistry but not for 10X processing.
    * `share.preprocess_tenx.barcode_dist`: When correcting 10X barcodes, this value defines the maximum number of mismatches allowed between the barcode and the whitelist `[default = 2]`.
    * `share.preprocess_tenx.threshold_pct_barcode_matching`: Defines the fraction of barcodes passing the barcode correction step required to consider the sample high-quality `[default = 0.6]`.
        >IMPORTANT: A threshold of `0.6` might be too high for your sample and usually is the main cause of failure when processing 10x samples.
    * It is possible to pass a custom barcode whitelist by setting the following parameter `share.whitelist`. In cases where RNA and ATAC use two different whitelists, such as `10x_multiome`, it is possible to set `share.whitelist_atac` and `share.whitelist_rna`.

6) [Resources](#resources)
    * It is not recommended to change the following parameters unless you need to increase the resources allocated for a certain task due to resource-related errors.

## Input files
* Input example for `shareseq`
    ```javascript
    {
        "share.read1_rna" : ["share-rna.R1.L1.fq.gz", "share-rna.R1.L2.fq.gz"],
        "share.read2_rna" : ["share-rna.R2.L1.fq.gz", "share-rna.R2.L2.fq.gz"],
        "share.read1_atac" : ["share-atac.R1.L1.fq.gz", "share-atac.R1.L2.fq.gz"],
        "share.read2_atac" : ["share-atac.R2.L1.fq.gz", "share-atac.R2.L2.fq.gz"],
    }
    ```
    * It is possible to run only the RNA or the ATAC portion of the experiment by leaving one of the two inputs empty. For example, to run only ATAC use the following:
    ```javascript
    {
        "share.read1_rna" : [],
        "share.read2_rna" : [],
        "share.read1_atac" : ["share-atac.R1.L1.fq.gz", "share-atac.R1.L2.fq.gz"],
        "share.read2_atac" : ["share-atac.R2.L1.fq.gz", "share-atac.R2.L2.fq.gz"],
    }
    ```
* Input example for `10x`. The pipeline expects in input FASTQs generated using the `cellranger mkfastq` command.
    ```javascript
    {
        "share.read1_rna" : ["10x-rna.R1.L1.fq.gz", "10x-rna.R1.L2.fq.gz"],
        "share.read2_rna" : ["10x-rna.R2.L1.fq.gz", "10x-rna.R2.L2.fq.gz"],
        "share.read1_atac" : ["10x-atac.R1.L1.fq.gz", "10x-atac.R1.L2.fq.gz"],
        "share.fastq_barcode" : ["10x-atac.R2.L1.fq.gz", "10x-atac.R2.L2.fq.gz"],
        "share.read2_atac" : ["10x-atac.R3.L1.fq.gz", "10x-atac.R3.L2.fq.gz"],
    }
    ```

## Resources

> **WARNING**: It is not recommended to change the following parameters unless you need to increase the resources allocated for a certain task due to resource-related errors.

For each task of the pipeline it is possible to fine-tune the resources used by the task using the following parameters: 

* `share.atac.[task_name]_cpus`
* `share.atac.[task_name]_memory_factor`
* `share.atac.[task_name]_disk_factor`
* `share.atac.[task_name]_docker_image`


For example, to change `align` runtimes parameters modify the following inputs:
* `share.atac.align_cpus`
* `share.atac.align_memory_factor`
* `share.atac.align_disk_factor`
* `share.atac.align_docker_image`
