version 1.0

struct Fastq {
	String pkrId
	String library
	String R1
	String sampleType
	String genome
	String notes
	Array[File] read1
	Array[File] read2
}

workflow wf_preprocess {
	input {
		# Preprocess inputs
		File bcl
		Boolean zipped = true
		Array[Int]? lanes
		File metaCsv
		String terra_project # set to none or make optional
		String workspace_name
		String dockerImage = "nchernia/share_task_preprocess:12"
	}

	String barcodeStructure = "14S10M28S10M28S9M8B"
	String sequencingCenter = "BI"
	String tar_flags = if zipped then 'xzf' else 'xf'
	String untarBcl =
		'gsutil -m -o GSUtil:parallel_thread_count=1' +
		' -o GSUtil:sliced_object_download_max_components=8' +
		' cp "~{bcl}" . && ' +
		'tar "~{tar_flags}" "~{basename(bcl)}" --exclude Images --exclude Thumbnail_Images' 
	
	String getSampleSheet =
		'gsutil -m -o GSUtil:parallel_thread_count=1' +
		' -o GSUtil:sliced_object_download_max_components=8' +
		' cp "~{bcl}" . && ' +
		'tar "~{tar_flags}" "~{basename(bcl)}" SampleSheet.csv'
	

	call BarcodeMap {
		input:
			metaCsv = metaCsv,
	}

	if (!defined(lanes)){
		call GetLanes { 
			input: 
				bcl = bcl, 
				untarBcl = getSampleSheet
		}
	}
	
	scatter (lane in select_first([lanes, GetLanes.lanes])) {
		call ExtractBarcodes {
			input:
				bcl = bcl,
				untarBcl = untarBcl,
				libraryBarcodes = BarcodeMap.out,
				barcodeStructure = barcodeStructure,
				lane = lane,
				dockerImage = dockerImage
		}

		call BasecallsToBams {
			input:
				bcl = bcl,
				untarBcl = untarBcl,
				barcodes = ExtractBarcodes.barcodes,
				libraryBarcodes = BarcodeMap.out,
				readStructure = ExtractBarcodes.readStructure,
				lane = lane,
				sequencingCenter = sequencingCenter,
				dockerImage = dockerImage
		}

		scatter(bam in BasecallsToBams.bams){
			# Convert unmapped, library-separated bams to fastqs
			# will assign cell barcode to read name 
			# assigns UMI for RNA to read name and adapter trims for ATAC
			call BamLookUp {
				input:
					bam = basename(bam),
					metaCsv = metaCsv,
			}

			call BamToFastq { 
				input: 
					bam=bam,
					pkrId = BamLookUp.pkrId,
					library = BamLookUp.library,
					sampleType = BamLookUp.sampleType, 
					genome = BamLookUp.genome,
					notes = BamLookUp.notes,
					R1barcodeSet = BamLookUp.R1barcodeSet, 
					dockerImage = dockerImage
			}

			call WriteTsvRow {
				input:
					fastq = BamToFastq.out
			}
		}
	}

	call QC {
		input:
			barcodeMetrics = ExtractBarcodes.barcodeMetrics
	}

	call GatherOutputs {
		input:
			rows = flatten(WriteTsvRow.row),
			name =  if zipped then basename(bcl, ".tar.gz") else basename(bcl, ".tar"),
			metaCsv = metaCsv, 
			dockerImage = dockerImage
	}

	call TerraUpsert {
		input:
			rna_tsv = GatherOutputs.rna_tsv,
			rna_no_tsv = GatherOutputs.rna_no_tsv,
			atac_tsv = GatherOutputs.atac_tsv,
			run_tsv = GatherOutputs.run_tsv,
			terra_project = terra_project,
			workspace_name = workspace_name, 
			dockerImage = dockerImage
	} 

	output {
		Array[String] percentMismatch = QC.percentMismatch
		Array[String] terraResponse = TerraUpsert.upsert_response
                Array[File] monitoringLogsExtract = ExtractBarcodes.monitoringLogExtract
                Array[File] monitoringLogsBasecalls = BasecallsToBam.monitoringLogBasecalls		
                # Array[Fastq] fastqs = flatten(BamToFastq.out)
		# Array[Array[Array[File]]] fastqs = BamToFastq.fastqs
	}
}

task BarcodeMap {
	input {
		File metaCsv
	}

	command <<<
		tail -n +6 ~{metaCsv} | cut -d, -f2 |  sed 's/ /\t/' > barcodes.tsv
	>>>

	output {
		Map[String, String] out = read_map("barcodes.tsv")
	}

	runtime {
		docker: "ubuntu:latest"
	}
}

task GetLanes {
	input {
		File bcl
		String untarBcl
	}

	parameter_meta {
		bcl: {
			localization_optional: true
		}
	}
	
	Float bclSize = size(bcl, 'G')
	Int diskSize = ceil(2.1 * bclSize)
	String diskType = if diskSize > 375 then "SSD" else "LOCAL"
	# Float memory = ceil(5.4 * bclSize + 147) * 0.25

	command <<<
		set -e

		~{untarBcl}
		tail -n+2 SampleSheet.csv | cut -d, -f2
	>>>

	output {
		Array[Int] lanes = read_lines(stdout())
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
		disks: "local-disk ~{diskSize} ~{diskType}"
		# memory: memory + 'G'
		# cpu: 14
	}
}

task ExtractBarcodes {
	input {
		# This function calls Picard to do library demultiplexing
		File bcl
		String untarBcl
		Map[String,String] libraryBarcodes
		String barcodeStructure 
		Int lane 
		String dockerImage
	}
	
	parameter_meta {
		bcl: {
			localization_optional: true
		}
	}

	File barcodesMap = write_map(libraryBarcodes)

	Int nBarcodes = 1
	String barcodeParamsFile = "barcode_params.tsv"
	String barcodeMetricsFile = "barcode_metrics.tsv"

	Float bclSize = size(bcl, 'G')

	Int diskSize = ceil(2.1 * bclSize)
	String diskType = if diskSize > 375 then "SSD" else "LOCAL"

	Float memory = ceil(0.2 * bclSize) * 1.25# an unusual increase from 0.25 x for black swan
	Int javaMemory = ceil((memory - 0.5) * 1000)

        String laneUntarBcl = untarBcl + ' RunInfo.xml RTAComplete.txt RunParameters.xml Data/Intensities/s.locs Data/Intensities/BaseCalls/L00~{lane}  && rm "~{basename(bcl)}"'
	command <<<
		set -e
		bash software/monitor_script.sh > monitoring.log &
		~{laneUntarBcl}

		# append terminating line feed
		sed -i -e '$a\' ~{barcodesMap}

		readLength=$(xmlstarlet sel -t -v "/RunInfo/Run/Reads/Read/@NumCycles" RunInfo.xml | head -n 1)T
		readStructure=${readLength}"~{barcodeStructure}"${readLength}
		echo ${readStructure} > readStructure.txt

		printf "barcode_name\tbarcode_sequence1" | tee "~{barcodeParamsFile}"
		while read -r params; do	
			name=$(echo "${params}" | cut -d$'\t' -f1)
			barcodes=$(echo "${params}" | cut -d$'\t' -f2-)
			printf "\n%s\t%s" "${name}" "${barcodes}" | tee -a "~{barcodeParamsFile}"
		done < "~{barcodesMap}"

		# Extract barcodes, write to metrics file
		java -Xmx~{javaMemory}m -jar /software/picard.jar ExtractIlluminaBarcodes \
			-BASECALLS_DIR "Data/Intensities/BaseCalls" \
			-TMP_DIR . \
			-OUTPUT_DIR . \
			-BARCODE_FILE "~{barcodeParamsFile}" \
			-METRICS_FILE "~{barcodeMetricsFile}" \
			-READ_STRUCTURE "${readStructure}" \
			-LANE "~{lane}" \
			-NUM_PROCESSORS 0 \
			-COMPRESSION_LEVEL 1 \
			-GZIP true
	>>>

	runtime {
		docker: dockerImage
		disks: "local-disk ~{diskSize} ~{diskType}"
		memory: memory + 'G'
		cpu: 4
	}

	output {
		String readStructure = read_string("readStructure.txt")
		File barcodeMetrics = barcodeMetricsFile
		File barcodes = write_lines(glob("*_barcode.txt.gz"))
		File monitoringLogExtract = "monitoring.log"
	}
}

task BasecallsToBams {
	input {
		# This function calls Picard to do library demultiplexing
		File bcl
		String untarBcl
		File barcodes
		Map[String,String] libraryBarcodes
		String readStructure 
		Int lane
		String sequencingCenter
		String dockerImage
	}

	parameter_meta {
		bcl: {
			localization_optional: true
		}
	}

	File barcodesMap = write_map(libraryBarcodes)
	String runIdFile = 'run_id.txt'
	String flowcellIdFile = 'flowcell_id.txt'
	String instrumentIdFile = 'instrument_id.txt'

	Float bclSize = size(bcl, 'G')

	Int diskSize = ceil(2.1 * bclSize)
	String diskType = if diskSize > 375 then "SSD" else "LOCAL"

	Float memory = ceil(1.4 * bclSize + 147) * 0.25
	Int javaMemory = ceil((memory - 0.5) * 1000)
        String laneUntarBcl = untarBcl + ' RunInfo.xml RTAComplete.txt RunParameters.xml Data/Intensities/s.locs Data/Intensities/BaseCalls/L00~{lane}  && rm "~{basename(bcl)}"'
	command <<<
		set -e
		bash software/monitor_script.sh > monitoring.log &
		~{laneUntarBcl}
		time gsutil -m cp -I . < "~{barcodes}"
		
		# append terminating line feed
		sed -i -e '$a\' ~{barcodesMap}

		# extract run parameters
		get_param () {
			param=$(xmlstarlet sel -t -v "/RunInfo/Run/$1" RunInfo.xml)
			echo "${param}" | tee "$2"
		}
		RUN_ID=$(get_param "@Number" "~{runIdFile}")
		FLOWCELL_ID=$(get_param "Flowcell" "~{flowcellIdFile}")
		INSTRUMENT_ID=$(get_param "Instrument" "~{instrumentIdFile}")

		# prepare library parameter files
		LIBRARY_PARAMS="library_params.tsv"
		printf "SAMPLE_ALIAS\tLIBRARY_NAME\tOUTPUT\tBARCODE_1\n" | tee "${LIBRARY_PARAMS}"
		while read -r params; do	
			name=$(echo "${params}" | cut -d$'\t' -f1)
			barcodes=$(echo "${params}" | cut -d$'\t' -f2-)
			printf "\n%s\t%s\t%s_L%d.bam\t%s" \
				"${name}" "${name}" "${name// /_}" "~{lane}" "${barcodes}" \
				| tee -a "${LIBRARY_PARAMS}"
		done < "~{barcodesMap}"

		# generate BAMs
		java -Xmx~{javaMemory}m -jar /software/picard.jar IlluminaBasecallsToSam \
			BASECALLS_DIR="Data/Intensities/BaseCalls" \
			BARCODES_DIR=. \
			TMP_DIR=. \
			LIBRARY_PARAMS="${LIBRARY_PARAMS}" \
			IGNORE_UNEXPECTED_BARCODES=true \
			INCLUDE_NON_PF_READS=false \
			READ_STRUCTURE="~{readStructure}" \
			LANE="~{lane}" \
			RUN_BARCODE="${INSTRUMENT_ID}:${RUN_ID}:${FLOWCELL_ID}" \
			SEQUENCING_CENTER="~{sequencingCenter}" \
			NUM_PROCESSORS=0 \
			MAX_RECORDS_IN_RAM=5000000 
	>>>

	runtime {
		docker: dockerImage
		disks: "local-disk ~{diskSize} ~{diskType}"
		memory: memory + 'G'
		cpu: 14
	}

	output {
		Array[File] bams = glob("*.bam")
                File monitoringLogBasecalls = "monitoring.log"
	}
}

task BamLookUp {
	# Find pkrId, sampleType, and barcodeSets from CSV
	# Rigid assumption about order of columns in CSV
	input {
		String bam
		File metaCsv
	}

	command <<<
		bucket="gs://broad-buenrostro-bcl-outputs/"
		file=~{bam}
		lib="${file%_*} "
		grep $lib ~{metaCsv} | cut -d, -f1 | sed 's/ /-/' > pkrId.txt
		echo ${file%_*} > library.txt
		barcode1=$(grep $lib ~{metaCsv} | cut -d, -f3)
		echo ${bucket}${barcode1}.txt > R1barcodeSet.txt
		grep $lib ~{metaCsv} | cut -d, -f4 > sampleType.txt
		grep $lib ~{metaCsv} | cut -d, -f5 > genome.txt
		grep $lib ~{metaCsv} | cut -d, -f6 > notes.txt
	>>>

	output {
		String pkrId = read_string("pkrId.txt")
		String library = read_string("library.txt")
		String R1barcodeSet = read_string("R1barcodeSet.txt")
		String sampleType = read_string("sampleType.txt")
		String genome = read_string("genome.txt")
		String notes = read_string("notes.txt")
	}

	runtime {
		docker: "ubuntu:latest"
	}
}

task BamToFastq {
	# Convert unmapped, library-separated bams to fastqs
	# will assign cell barcode to read name 
	# assigns UMI for RNA to read name and adapter trims for ATAC

	# Defaults to file R1.txt in the src/python directory if no round barcodes given
	input {
		File bam
		String pkrId
		String library
		String sampleType
		String genome
		String notes
		# Array[Array[String]] R1barcodeSet
		# Array[Array[String]]? R2barcodes
		# Array[Array[String]]? R3barcodes
		File R1barcodeSet
		File? R2barcodes
		File? R3barcodes
		String dockerImage
	}

	String prefix = basename(bam, ".bam")
	
	Float bamSize = size(bam, 'G')

	Int diskSize = ceil(bamSize + 5)
	String diskType = if diskSize > 375 then "SSD" else "LOCAL"

	Float memory = ceil(1.5 * bamSize + 1) * 2


	# Workaround since write_tsv does not take type "?", must be defined
	# Array[Array[String]] R2_if_defined = select_first([R2barcodes, []])
	# Array[Array[String]] R3_if_defined = select_first([R3barcodes, []])

	# Use round 1 default barcode set in rounds 2 and 3 if not sent in
	File R1file = R1barcodeSet #write_tsv(R1barcodeSet)
	File R2file = if defined(R2barcodes)
					# then write_tsv(R2_if_defined) else R1file
					then R2barcodes else R1barcodeSet
	File R3file = if defined(R3barcodes)
					# then write_tsv(R3_if_defined) else R1file
					then R3barcodes else R1barcodeSet


	command <<<
		python3 /software/bam_fastq.py ~{bam} ~{R1file} ~{R2file} ~{R3file} ~{pkrId} -p ~{prefix} -s ~{sampleType}

		gzip *.fastq
	>>>

	output {
		Fastq out = object {
			pkrId: pkrId,
			library: library,
			R1: R1barcodeSet,
			sampleType: sampleType,
			genome: genome,
			notes: notes,
			read1: glob("*R1.fastq.gz"),
			read2: glob("*R2.fastq.gz")
		}
		# Array[File] fastqs = glob("*.fastq")
	}
	runtime {
		docker: dockerImage
		disks: "local-disk ~{diskSize} ~{diskType}"
		memory: memory + 'G'
	}
}

task QC {
	input {
		Array[File] barcodeMetrics
	}

	Int total = length(barcodeMetrics)
	
	command <<<
		ARRAY=(~{sep=" " barcodeMetrics}) # Load array into bash variable
		for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;)
		do
			awk '$1=="NNNNNNNN"' ${ARRAY[$c]} | cut -f11
		done
	>>>
	
	output {
		Array[String] percentMismatch = read_lines(stdout())
	}
	
	runtime {
		docker: "ubuntu:latest"
	}
}


task WriteTsvRow {
	input {
		Fastq fastq
	}

	Array[String] read1 = fastq.read1
	Array[String] read2 = fastq.read2

	command <<<
		# echo -e "Library\tPKR\tR1_subset\tType\tfastq_R1\tfastq_R2\tGenome\tNotes" > fastq.tsv
		echo -e "~{fastq.library}\t~{fastq.pkrId}\t~{fastq.R1}\t~{fastq.sampleType}\t~{sep=',' read1}\t~{sep=',' read2}\t~{fastq.genome}\t~{fastq.notes}" > row.tsv
	>>>

	output {
		File row = 'row.tsv'
	}

	runtime {
		docker: "ubuntu:latest"
	}
}

task GatherOutputs {
	input {
		Array[File] rows
		String name
		File metaCsv
		String dockerImage
	}

	command <<<
		echo -e "Library\tPKR\tR1_subset\tType\tfastq_R1\tfastq_R2\tGenome\tNotes" > fastq.tsv
		cat ~{sep=' ' rows} >> fastq.tsv

		python3 /software/write_terra_tables.py --input 'fastq.tsv' --name ~{name} --meta ~{metaCsv}
	>>>

	runtime {
		docker: dockerImage
	}
	output {
		File rna_tsv = "rna.tsv"
		File rna_no_tsv = "rna_no.tsv"
		File atac_tsv = "atac.tsv"
		File run_tsv = "run.tsv"
	}
}

task TerraUpsert {
	input {
		File rna_tsv
		File rna_no_tsv
		File atac_tsv
		File run_tsv
		String terra_project
		String workspace_name
		String dockerImage
	}
	
	command <<<
		set -e
		python3 /software/flexible_import_entities_standard.py \
			-t "~{rna_tsv}" \
			-p "~{terra_project}" \
			-w "~{workspace_name}"
		
		python3 /software/flexible_import_entities_standard.py \
			-t "~{rna_no_tsv}" \
			-p "~{terra_project}" \
			-w "~{workspace_name}"

		python3 /software/flexible_import_entities_standard.py \
			-t "~{atac_tsv}" \
			-p "~{terra_project}" \
			-w "~{workspace_name}"

		python3 /software/flexible_import_entities_standard.py \
			-t "~{run_tsv}" \
			-p "~{terra_project}" \
			-w "~{workspace_name}"
	>>>
	
	runtime {
		docker: dockerImage
		memory: "2 GB"
		cpu: 1
	}
	
	output {
		Array[String] upsert_response = read_lines(stdout())
	}
}
