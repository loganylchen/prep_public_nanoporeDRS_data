rule slow5tools_f2s:
    input:
        raw_fast5_dir = "data/tmp.{sample}.nanoporeDRS",
        tag = "data/tmp.{sample}.nanoporeDRS.tag",
        check_tag="data/tmp.{sample}.nanoporeDRS.check.tag"
    output:
        raw_slow5_dir = temp("data/tmp.{sample}.slow5"),
        tag=temp("data/tmp.{sample}.slow5.tag"),
    threads: config["threads"]["slow5tools"]
    conda:
        "../envs/slow5tools.yaml"
    log:
        "logs/slow5tools_f2s/{sample}.log"
    benchmark:
        "logs/slow5tools_f2s/{sample}.benchmark"
    shell:
        'slow5tools f2s -p {threads} -a -d {input.raw_fast5_dir}  {output.raw_slow5_dir} 2>{log} && '
        'touch {output.tag}'

rule slow5tools_merge:
    input:
        raw_slow5_dir="data/tmp.{sample}.slow5",
        tag="data/tmp.{sample}.slow5.tag",
    output:
        blow5="data/{sample}/blow5/nanopore.blow5"
    threads: config["threads"]["slow5tools"]
    conda:
        "../envs/slow5tools.yaml"
    log:
        "logs/slow5tools_merge/{sample}.log"
    benchmark:
        "logs/slow5tools_merge/{sample}.benchmark"
    shell:
        "slow5tools merge -t {threads} -o {output.blow5} {input.raw_slow5_dir} 2>{log}"


