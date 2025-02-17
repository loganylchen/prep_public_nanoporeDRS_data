rule basecalling:
    input:
        raw_fast5_dir="data/tmp.{sample}.nanoporeDRS",
        tag="data/tmp.{sample}.nanoporeDRS.tag",
        check_tag="data/tmp.{sample}.nanoporeDRS.check.tag",
    output:
        raw_fastq_dir=temp(directory("data/tmp.{sample}.fastq")),
        tag=temp("data/tmp.{sample}.fastq.tag"),
    params:
        guppy=config['guppy']['path'],
        ext_param=config['guppy']['param']
    threads: config['threads']['guppy']
    log:
        log="logs/basecalling/{sample}.log",
        err="logs/basecalling/{sample}.err",
    benchmark:
        "benchmarks/basecalling/{sample}.txt"
    shell:
        '{params.guppy} -i {input.raw_fast5_dir} "
        '-s {output.raw_fastq_dir} '
        '{params.ext_param} 2>>{log.err} && touch {output.tag}'


rule merging:
    input:
        raw_fastq_dir="data/tmp.{sample}.fastq",
        tag="data/tmp.{sample}.fastq.tag"
    output:
        fastq="{project}/data/{sample}/fastq/pass.fq.gz"
    threads: 1
    log:
        "logs/merging/{sample}_{project}.log",
    benchmark:
        "benchmarks/basecalling/{sample}_{project}.txt"
    shell:
        'zcat {input.raw_fastq_dir}/pass/*fastq.gz | gzip -c > {output.fastq}'


