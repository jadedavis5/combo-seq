// Module information
name = "seqkit"
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
if (module.memory == -1) {
    throw new Exception(
        "Module ${name} does not support dynamic memory assignment, \
        please set to 0 or a positive integer")
} else if (module.memory == 0) {
    memory = 1
    module.memory = "${module.memory}G"
    binding.variables.remove 'memory'
}
else {
    module.memory = "${module.memory}G"
}


// Module settings
logfile = ""


process length {
    tag "${name}-${seq}"
    cpus module.cores
    memory module.memory
    executor "local"
    // time module.time
    // queue module.queue

    // publishDir "${params.data.out}/${params.data.time}",
    //     mode: "copy",
    //     overwrite: true,
    //     pattern: "time-${id}-${name}.txt"
    // publishDir "${params.data.logs}/${params.data.star}",
    //     mode: "copy",
    //     overwrite: true,
    //     pattern: "${logfile}",
    // publishDir "${params.data.logs}/${params.data.star}",
    //     mode: "copy",
    //     overwrite: true,
    //     pattern: "${logfile}",
    //     saveAs: {"${name}-index-${id}-${it}"}
    //     saveAs: {"${name}-index-${id}-${it}"}

    input:
    path(seq)

    output:
    env(seq_length), emit: "seq_length"

    shell:
    '''
    seq_length=$(!{bin} stats -T "!{seq}" -j !{task.cpus} |
        cut -d $'\n' -f 2 |
        cut -d $'\t' -f 5 )
    '''
}
