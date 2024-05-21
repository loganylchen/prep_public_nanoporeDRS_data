rule fetch_data:
    output:
        compressed_data = temp("data/tmp.{sample}.nanoporeDRS.tar.gz"),
        tag = temp("data/tmp.{sample}.nanoporeDRS.tar.gz.tag")
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
        compressed_data = "data/tmp.{sample}.nanoporeDRS.tar.gz",
        tag= "data/tmp.{sample}.nanoporeDRS.tar.gz.tag"
    output:
        extracted_data = temp("data/tmp.{sample}.nanoporeDRS"),
        tag = temp("data/tmp.{sample}.nanoporeDRS.tag")
    log:
        "logs/extract_data/{sample}.log"
    benchmark:
        "benchmarks/extract_data/{sample}.txt"
    shell:
        "mkdir -p {output.extracted_data} && "
        "tar -xzf {input.compressed_data} -C {output.extracted_data} 2>{log} && touch {output.tag}"


