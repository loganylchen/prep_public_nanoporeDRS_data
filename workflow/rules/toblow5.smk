rule slow5tools_f2s:
    input:
        raw_fast5_dir = "data/tmp.{sample}.nanoporeDRS",
        tag = "data/tmp.{sample}.nanoporeDRS.tag",
        check_tag="data/tmp.{sample}.nanoporeDRS.check.tag"
    output:
        raw_slow5_dir = temp(directory("data/tmp.{sample}.slow5")),
        tag=temp("data/tmp.{sample}.slow5.tag"),
    threads: config["threads"]["slow5tools"]
    container:
        "btrspg/slow5tools:latest"
    log:
        log="logs/slow5tools_f2s/{sample}.log",
        err="logs/slow5tools_f2s/{sample}.err"
    benchmark:
        "benchmarks/slow5tools_f2s/{sample}.benchmark"
    shell:
        'slow5tools f2s -p {threads} -a -d {output.raw_slow5_dir} {input.raw_fast5_dir}  1>{log.log} 2>{log.err} && '
        'touch {output.tag}'

rule slow5tools_merge:
    input:
        raw_slow5_dir="data/tmp.{sample}.slow5",
        tag="data/tmp.{sample}.slow5.tag",
    output:
        blow5="data/{sample}/blow5/nanopore.blow5"
    threads: config["threads"]["slow5tools"]
    container:
        "btrspg/slow5tools:latest"
    log:
        log="logs/slow5tools_merge/{sample}.log",
        err="logs/slow5tools_merge/{sample}.err"
    benchmark:
        "benchmarks/slow5tools_merge/{sample}.benchmark"
    shell:
        "slow5tools merge -t {threads} -o {output.blow5} {input.raw_slow5_dir} 1>{log.log} 2>{log.err}"


