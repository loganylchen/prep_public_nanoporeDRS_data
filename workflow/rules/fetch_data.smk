rule fetch_data:
    output:
        compressed_data=temp("data/tmp.{sample}.nanoporeDRS.data"),
        tag=temp("data/tmp.{sample}.nanoporeDRS.download.tag")
    params:
        url = get_data_url
    threads: config["threads"]["fetch_data"]
    conda:
        "../envs/fetch_data.yaml"
    log:
        "logs/fetch_data/{sample}.log"
    benchmark:
        "benchmarks/fetch_data/{sample}.txt"
    shell:
        'lftp -c "pget -n {threads} {params.url} -o {output.compressed_data}" >{log} && touch {output.tag}'

rule extract_data:
    input:
        compressed_data="data/tmp.{sample}.nanoporeDRS.data",
        tag="data/tmp.{sample}.nanoporeDRS.download.tag"
    output:
        extracted_data = temp(directory("data/tmp.{sample}.nanoporeDRS")),
        tag = temp("data/tmp.{sample}.nanoporeDRS.tag")
    params:
        command = get_uncompress_command
    log:
        "logs/extract_data/{sample}.log"
    benchmark:
        "benchmarks/extract_data/{sample}.txt"
    shell:
        "mkdir -p {output.extracted_data} && "
        "{params.command} {input.compressed_data} -C {output.extracted_data} 2>{log} && touch {output.tag}"

rule check_data:
    input:
        extracted_data = "data/tmp.{sample}.nanoporeDRS",
        tag = "data/tmp.{sample}.nanoporeDRS.tag"
    output:
        tag= temp("data/tmp.{sample}.nanoporeDRS.check.tag")
    log:
        "logs/check_data/{sample}.log"
    benchmark:
        "benchmarks/check_data/{sample}.txt"
    script:
        "../scripts/check_data.py"


