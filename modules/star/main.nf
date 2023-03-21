// Module information
name = "star"
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
// Get read length
overhang = params.data.reads_length.toInteger() -
    params.fastp.trimf.toInteger() -
    params.fastp.trimt.toInteger() - 1
logfile = "Log.out"


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
    tuple val(id), path(reads), path(genome), val(genome_length), path(gtf), val(genome_name), path("*")

    output:
    stdout
    tuple val(id), emit: "star"
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true

    shell:
    '''
    flags="!{module.flags}"
    if [[ "!{params.workflow.compress}" == "True" ]]; then
        flags+=" --readFilesCommand gunzip -c"
    fi

    # Disable shared memory when using Slurm, otherwise enable
    if [[ "!{task.executor}" == "slurm" ]]; then
        flags+=" --genomeLoad NoSharedMemory"
        shmem=false
    else
        # Disable for two-pass mapping and on-the-fly junction insertion
        # flags+=" --genomeLoad LoadAndKeep"
        shmem=false
    fi

    # Calculate STAR suffix array index
    x="!{genome_length}"
    x=$(echo "scale = 10; (l(${x})/l(2) - 1)/2 - 1" | "!{params.tools.bc}" -l)
    x=$(echo "scale = 0; ${x}/1" | "!{params.tools.bc}")
    if [[ x -gt 14 ]]; then x=14; fi

    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} ${flags} \
        --readFilesIn "!{reads[0]}" "!{reads[1]}" \
        --genomeFastaFiles "!{genome}" \
        --genomeDir "!{genome_name}" \
        --quantMode TranscriptomeSAM GeneCounts \
        --sjdbGTFfile "!{gtf}" \
        --sjdbOverhang "!{overhang}" \
        --sjdbInsertSave "All" \
        --genomeSAindexNbases ${x} \
        --runThreadN "!{task.cpus}" \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN "!{task.cpus}" \
        --outReadsUnmapped Fastx

    # Remove shared memory
    if [[ ${shmem} == true ]]; then
        !{bin} --genomeLoad Remove
    fi

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}

process index {
    tag "${name}-index-${id}"
    cpus module.cores
    memory module.memory
    time module.time
    queue module.queue

    // pattern: "{Genome,SA,SAindex,chr*,exon*,gene*,genome*,sjdb*,transcript*}"
    publishDir "${params.data.out}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "${id}"
    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"
    publishDir "${params.data.logs}/${params.data.star}",
        mode: "copy",
        overwrite: true,
        pattern: "${logfile}",
        saveAs: {"${name}-index-${id}-${it}"}

    input:
    tuple val(id), path(genome), val(genome_length), path(gtf)

    output:
    // path("{Genome,SA,SAindex,chr*,exon*,gene*,genome*,sjdb*,transcript*}")
    tuple val(id),
        path("${id}"),
        emit: "index"
    path("${logfile}"), emit: "logs", optional: true

    shell:
    '''
    flags="!{module.flags}"

    # Calculate STAR suffix array index
    x="!{genome_length}"
    x=$(echo "scale = 1; (l(${x})/l(2) - 1)/2 - 1" | "!{params.tools.bc}" -l)
    x=$(echo "scale = 0; ${x}/1" | "!{params.tools.bc}")
    if [[ x -gt 14 ]]; then x=14; fi

    # Calculate STAR memory limits
    # Set to 16GB + 2 * genome_length
    y=$(echo "scale = 0; (8 * 1024^3) + (!{genome_length} * 2.5)/1" | "!{params.tools.bc}")

    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} "${flags}" \
        --genomeFastaFiles "!{genome}" \
        --genomeDir "!{id}" \
        --quantMode GeneCounts \
        --sjdbGTFfile "!{gtf}" \
        --sjdbOverhang "!{overhang}" \
        --sjdbInsertSave All \
        --genomeSAindexNbases ${x} \
        --runThreadN "!{task.cpus}" \
        --runMode genomeGenerate \
        --limitGenomeGenerateRAM ${y}

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
