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
    if (module.direct == "True") {
        bin = "${container}"
    } else {
        bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
    }
} else {
    bin = "${module.bin}"
}

// Module settings
flags = ""
if (module.autodetect == "True") {
    flags += " --detect_adapter_for_pe"
}
flags = flags.strip()


process module {
    tag "${name}-${id}"
    cpus module.cpus
    memory module.memory
    time module.time
    queue module.queue

    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: {"${tfile}"}

    input:
    tuple val(id), path(reads), path(adapters)

    output:
    tuple val(id), path("sub-${id}*")
    tuple env(tfile), path("${tfile}"), emit: "tfile"

    shell:
    '''
    flags=""
    if [[ "!{module.autodetect}" == "True" ]]; then
        flags+="${flags} !{flags}"
    else
        flags+="${flags} --adapter_fasta !{adapters[0]}"
    fi

    if [[ "!{params.workflow.compress}" == "True" ]] && [[ "!{module.compress}" == "True" ]]; then
        flags+="${flags} -z !{module.compress_level}"
    fi

    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} ${flags} \
        --thread !{task.cpus} \
        -i "!{reads[0]}" -I "!{reads[1]}" \
        -o "trim-!{reads[0]}" -O "trim-!{reads[1]}" \
        -q "!{module.quality}"

    # Process information
    echo " 	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
