rule fetch_data:
    output:
        unpack(get_compressed_data)
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
        unpack(get_compressed_data_notemp),
    output:
        extracted_data = temp(directory("data/tmp.{sample}.nanoporeDRS")),
        tag = temp("data/tmp.{sample}.nanoporeDRS.tag")
    params:
        unpack(get_uncompress_command)
    log:
        "logs/extract_data/{sample}.log"
    benchmark:
        "benchmarks/extract_data/{sample}.txt"
    shell:
        "mkdir -p {output.extracted_data} && "
        "{params.uncompress} {input.compressed_data} -C {output.extracted_data} 2>{log} && touch {output.tag}"

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


