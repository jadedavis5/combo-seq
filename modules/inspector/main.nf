// Module information
name = "inspector"
// module = params[name]


// Module settings


process module {
    executor "local"
    // cpus module.cores
    // memory module.memory
    // time module.time
    // queue module.queue

    input:
    tuple val(id), path("*")

    output:
    stdout

    shell:
    '''
    find .
    '''
}
