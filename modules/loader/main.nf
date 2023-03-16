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
    tuple val(query), val(id), path(reads)

    output:
    tuple val(id), path("trim-${id}*"), emit: "reads"

    shell:
    '''
    find "!{params.data.path}/!{query}" \
        -iname "*trim-!{id}*" \
        -exec bash -c "ln -sr {} ." ';'
    '''
}
