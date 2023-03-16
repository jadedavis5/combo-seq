// Module information
name = "extract"
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
        bin = "${container} ${module.bin}"
    } else {
        bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
    }
} else {
    bin = "${module.bin}"
}


// Set memory
if (module.memory == "-1") {
    memory = 1 + Math.floor(module.cores.toInteger() * 0.25)
    module.memory = "${memory}G"
    binding.variables.remove 'memory'
} else if (module.memory == "0") {
    memory = 1
    module.memory = "${module.memory}G"
    binding.variables.remove 'memory'
}
else {
    module.memory = "${module.memory}G"
}


// Module settings


process module {
    tag "${name}-${id}"
    cpus module.cores
    memory module.memory
    time module.time
    queue module.queue

    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}*.fastq"), emit: "reads"
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true

    shell:
    '''
    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{params.parallel.bin} '!{bin} -p !{task.cpus} -dc {} > {.}' ::: *.gz

    # Process information
    echo " 	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
