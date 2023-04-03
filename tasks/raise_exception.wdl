version 1.0

task exception {
    input {
        Boolean fail
    }
    command <<<
        >&2 echo 'Hello!!!! }{}"'
        if [ "~{fail}" == "true" ]; then
            echo '{"wdl_error_message": "this is the end, my only friend, the end", "meaning": 420}' > wdl_failure_message.json
            exit 42
        fi
    >>>

   runtime {
     failureMessageFile: "wdl_failure_message.json"
   }
}
