// Module information
name = "template"
module = params[name]

// Enable containers
if (module.container == true) {
    if (params.container.monolithic == true) {
        container = params.container.path
    } else {
        container = module.path
    }
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
        echo

    # Process information
    echo " 	Allocated resources" >> "$TFILE"
    echo "	CPUs: !{task.cpus}" >> "$TFILE"
    echo "	Memory: !{task.memory}" >> "$TFILE"
    echo "	Time: !{task.time}" >> "$TFILE"
    echo "	Queue: !{task.queue}" >> "$TFILE"
    '''
}
