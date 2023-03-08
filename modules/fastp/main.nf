// Module information
name = "fastp"
module = params[name]

// Enable containers
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

process module {
    tag "${name}-${id}"
    cpus module.cpus
    memory module.memory
    time module.time
    queue module.queue

    input:
    tuple val(id)

    output:
    tuple val(id)

    shell:
    '''
    TFILE="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "$TFILE" \
        !{bin} --help

    # Process information
    echo " 	Allocated resources" >> "$TFILE"
    echo "	CPUs: !{task.cpus}" >> "$TFILE"
    echo "	Memory: !{task.memory}" >> "$TFILE"
    echo "	Time: !{task.time}" >> "$TFILE"
    echo "	Queue: !{task.queue}" >> "$TFILE"
    '''
}
