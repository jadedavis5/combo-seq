// Module information
name = "inspector"
// module = params[name]


// Module settings


process module {
    executor "local"
    // cpus module.cpus
    // memory module.memory
    // time module.time
    // queue module.queue

    input:
    tuple val(x)

    output:
    stdout

    shell:
    '''
    echo "!{x}"
    '''
}
