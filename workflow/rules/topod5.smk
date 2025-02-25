rule pod5_f2p:
    input:
        raw_fast5_dir="data/tmp.{sample}.nanoporeDRS",
        tag="data/tmp.{sample}.nanoporeDRS.tag",
        check_tag="data/tmp.{sample}.nanoporeDRS.check.tag",
    output:
        pod5="{project}/data/{sample}/pod5/nanopore.pod5",
        
    threads: config["threads"]["pod5"]
    conda:
        "../envs/pod5.yaml"
    priority: 9
    log:
        log="logs/pod5_f2p/{sample}.log",
        err="logs/pod5_f2p/{sample}.err",
    benchmark:
        "benchmarks/pod5_f2p/{sample}.benchmark"
    shell:
        "pod5 convert fast5 -t {threads} -f -r -o {output.pod5} {input.raw_fast5_dir}  1>{log.log} 2>{log.err} "


# rule pod5_merge:
#     input:
#         raw_pod5_dir="data/tmp.{sample}.pod5" ,
#         tag="data/tmp.{sample}.pod5.tag",
#     output:
#         pod5="{project}/data/{sample}/pod5/nanopore.pod5",
#     threads: config["threads"]["pod5"]
#     conda:
#         "../envs/pod5.yaml"
#     priority: 10
#     log:
#         log="logs/pod5_merge/{sample}_{project}.log",
#         err="logs/pod5_merge/{sample}_{project}.err",
#     benchmark:
#         "benchmarks/pod5_merge/{sample}_{project}.benchmark.txt"
#     shell:
#         "pod5 merge -t {threads} -o {output.pod5} -f {input.raw_pod5_dir} 1>{log.log} 2>{log.err}"
