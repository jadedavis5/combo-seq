// Module information
name = "fastp"
module = params[name]

// Enable containers and set binary
if (module.container == "True") {
    // Check if monolithic containers are enabled, and set container
    if (params.container.monolithic == "True") {
        container = params.container.path
    } else {
        container = module.path
    }

    // Set executable
    bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
} else {
    bin = "${module.bin}"
}

println(module.autodetect)
if (module.autodetect == "True") {
    flags += " --detect_adapter_for_pe"
} else {

}
flags = flags.strip()

process module {
    tag "${name}-${id}"
    cpus module.cpus
    memory module.memory
    time module.time
    queue module.queue

    input:
    tuple val(id), path(reads), path(adapters)

    output:
    tuple val(id)

    shell:
    '''
    if [[ "!{module.autodetect}" != "True" ]]; then
    flags+=" --adapter_fasta !{adapters[0]}"
    fi

    TFILE="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "$TFILE" \
        !{bin} !{flags} \
        -i "!{reads[0]}" -I "!{reads[1]}" \
        -o "trim-!{reads[0]}" -O "trim-!{reads[1]}" \
        -q "!{module.quality}"

    # Process information
    echo " 	Allocated resources" >> "$TFILE"
    echo "	CPUs: !{task.cpus}" >> "$TFILE"
    echo "	Memory: !{task.memory}" >> "$TFILE"
    echo "	Time: !{task.time}" >> "$TFILE"
    echo "	Queue: !{task.queue}" >> "$TFILE"
    '''
}
