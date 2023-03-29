// Module information
name = "loader"


process loader1 {
    tag "${name}"
    cpus 1
    memory "1G"
    executor "local"
    // time module.time
    // queue module.queue

    input:
    tuple val(id), val(query)

    output:
    tuple val(id), path("${id}*")

    shell:
    '''
    find "!{params.data.path}/!{query}" \
        -iname "*!{id}*" \
        -exec bash -c "ln -sr {} ." ';'
    '''
}

process loader2 {
    tag "${name}"
    cpus 1
    memory "1G"
    executor "local"
    // time module.time
    // queue module.queue

    input:
    tuple val(id), val(path), val(query)

    output:
    tuple val(id), path("${query}")

    shell:
    '''
    find "!{params.data.path}/!{path}/!{id}" \
        -iname "!{query}" \
        -exec bash -c "ln -sr {} ." ';'
    '''
}
