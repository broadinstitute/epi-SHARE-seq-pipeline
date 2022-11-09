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
    	docker : 'encodedcc/chip-seq-pipeline:v2.2.1'
    }
}
