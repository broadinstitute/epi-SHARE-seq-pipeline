version 1.0

workflow wf_preprocess {
  input {
    # Preprocess inputs
    File bcl
    Map[String,String] libraryBarcodes
    String readStructure = "50T14S10M28S10M28S9M8B50T"
    String sequencingCenter = "BI"
    String dockerImage = "nchernia/share_task_preprocess:2"
  }

  String untarBcl =
    'gsutil -m -o GSUtil:parallel_thread_count=1' +
    ' -o GSUtil:sliced_object_download_max_components=8' +
    ' cp "~{bcl}" . && ' +
    'tar xf "~{basename(bcl)}" --exclude Images --exclude Thumbnail_Images && ' +
    'rm "~{basename(bcl)}"'

  call GetLanes { 
    input: bcl = bcl, untarBcl = untarBcl
  }

  scatter (lane in GetLanes.lanes) {
   call ExtractBarcodes {
      input:
        bcl = bcl,
        untarBcl = untarBcl,
        libraryBarcodes = libraryBarcodes,
        readStructure = readStructure,
        lane = lane,
        dockerImage = dockerImage
    }

    call BasecallsToBams {
        input:
          bcl = bcl,
          untarBcl = untarBcl,
          barcodes = ExtractBarcodes.barcodes,
          libraryBarcodes = libraryBarcodes,
          readStructure = readStructure,
          lane = lane,
          sequencingCenter = sequencingCenter,
          dockerImage = dockerImage
      }
  }
  output {
    Array[Array[File]] bams = BasecallsToBams.bams
	Array[File] barcodeMetrics = ExtractBarcodes.barcodeMetrics
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

  Int diskSize = ceil(1.9 * bclSize + 5)
  String diskType = if diskSize > 375 then "SSD" else "LOCAL"
  Float memory = ceil(5.4 * bclSize + 147) * 0.25
  command {
    set -e

    ~{untarBcl}
    tail -n+2 SampleSheet.csv | cut -d, -f2
  }
  output {
    Array[Int] lanes = read_lines(stdout())
  }
  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
    disks: "local-disk ~{diskSize} ~{diskType}"
    memory: memory + 'G'
    cpu: 14
  }
}

task ExtractBarcodes {
   input {
    # This function calls Picard to do library demultiplexing
     File bcl
     String untarBcl
     Map[String,String] libraryBarcodes
     String readStructure 
     Int lane 
     String dockerImage
  }
  parameter_meta {
    bcl: {
      localization_optional: true
    }
  }
  String barcodeParamsFile = "barcode_params.tsv"
  String barcodeMetricsFile = "barcode_metrics.tsv"
  
  
  Float bclSize = size(bcl, 'G')

  Int diskSize = ceil(1.9 * bclSize + 5)
  String diskType = if diskSize > 375 then "SSD" else "LOCAL"

  Float memory = ceil(0.8 * bclSize) * 1.25  # an unusual increase from 0.25 x for black swan
  Int javaMemory = ceil((memory - 0.5) * 1000)
  Int nBarcodes = 1
  command <<<
    ~{untarBcl}
 
    printf "barcode_name\tbarcode_sequence1" | tee "~{barcodeParamsFile}"
    while read -r params; do	
      name=$(echo "${params}" | cut -d$'\t' -f1)
      barcodes=$(echo "${params}" | cut -d$'\t' -f2-)
      printf "\n%s\t%s" "${name}" "${barcodes}" | tee -a "~{barcodeParamsFile}"
    done < "~{write_map(libraryBarcodes)}"
    # Extract barcodes, write to metrics file
    java -Xmx~{javaMemory}m -jar /software/picard.jar ExtractIlluminaBarcodes \
      -BASECALLS_DIR "Data/Intensities/BaseCalls" \
      -TMP_DIR . \
      -OUTPUT_DIR . \
      -BARCODE_FILE "~{barcodeParamsFile}" \
      -METRICS_FILE "~{barcodeMetricsFile}" \
      -READ_STRUCTURE "~{readStructure}" \
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
    File barcodeMetrics = barcodeMetricsFile
    File barcodes = write_lines(glob("*_barcode.txt.gz"))
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
#    Array[Array[String]] round1BarcodeSet
  
  parameter_meta {
    bcl: {
      localization_optional: true
    }
  }
  
  String runIdFile = 'run_id.txt'
  String flowcellIdFile = 'flowcell_id.txt'
  String instrumentIdFile = 'instrument_id.txt'

  Float bclSize = size(bcl, 'G')

  Int diskSize = ceil(1.9 * bclSize + 5)
  String diskType = if diskSize > 375 then "SSD" else "LOCAL"

  Float memory = ceil(5.4 * bclSize + 147) * 0.25
  Int javaMemory = ceil((memory - 0.5) * 1000)
    

  command <<<
    set -e
    
    ~{untarBcl}
    time gsutil -m cp -I . < "~{barcodes}"
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
    done < "~{write_map(libraryBarcodes)}"
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
  }
}

#task gather_outputs {
#    input {
#        Array[String] read_paths
#        File yaml
#        String new_table_name       
#    }
#
#    command <<<
#        python3 /software/gather_outputs.py -n ~{new_table_name} -y ~{yaml} -p ~{sep="," read_paths}
#    >>>
#
#    runtime {
#        docker: "polumechanos/share-seq-undetermined"
#    }
#    output {
#        File load_tsv = "results_load_table.tsv"
#    }
#}
#task entities_batch_upsert {
#  input {
#    File      tsv_file
#    String    terra_project
#    String    workspace_name
#  }
#    String    docker="polumechanos/update_terra_table"
#  command {
#    set -e
#    python3 /software/batch_upsert_entities_standard.py \
#        -t "~{tsv_file}" \
#        -p "~{terra_project}" \
#        -w "~{workspace_name}"
#  }
#  runtime {
#    docker : docker
#    memory: "2 GB"
#    cpu: 1
#  }
#  output {
#    String upsert_entities_response = stdout()
#  }
#}
