// Module information
name = "loader"


process module {
    tag "${name}"
    cpus 1
    memory "1G"
    executor "local"
    // time module.time
    // queue module.queue

    input:
    tuple val(query), val(id)

    output:
    tuple val(id), path("${id}*")

    shell:
    '''
    find "!{params.data.path}/!{query}" \
        -iname "*!{id}*" \
        -exec bash -c "ln -sr {} ." ';'
    '''
}
