rule guppy_basecalling:
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
    priority: 9
    log:
        log="logs/guppy_basecalling/{sample}.log",
        err="logs/guppy_basecalling/{sample}.err",
    benchmark:
        "benchmarks/guppy_basecalling/{sample}.txt"
    shell:
        '{params.guppy} -i {input.raw_fast5_dir} '
        '-s {output.raw_fastq_dir} '
        '{params.ext_param} 2>>{log.err} && touch {output.tag}'


rule fastq_merge:
    input:
        raw_fastq_dir="data/tmp.{sample}.fastq" ,
        tag="data/tmp.{sample}.fastq.tag"
    output:
        fastq="{project}/data/{sample}/fastq/pass.fq.gz"
    threads: 1
    priority: 10
    log:
        log="logs/fastq_merge/{sample}_{project}.log",
    benchmark:
        "benchmarks/fastq_merge/{sample}_{project}.txt"
    shell:
        'zcat {input.raw_fastq_dir}/pass/*fastq.gz | gzip -c > {output.fastq}'


rule dorado_basecalling:
    input:
        pod5="{project}/data/{sample}/pod5/nanopore.pod5",
        model=f"tmp/models/{config['dorado']['model_name']}"
    output:
        raw_bam_dir=temp(directory("data/tmp.{sample}.bam")),
        tag=temp("data/tmp.{sample}.bam.tag"),
    params:
        dorado=config['dorado']['path'],
        model=config['dorado']['model_name'],
        mm2_opts=config['dorado']['mm2_opts'],
        reference=config['dorado']['reference'],
        ext_param=config['dorado']['param']
    threads: config['threads']['dorado']
    priority: 9
    log:
        log="logs/dorado_basecalling/{sample}.log",
        err="logs/dorado_basecalling/{sample}.err",
    benchmark:
        "benchmarks/dorado_basecalling/{sample}.txt"
    shell:
        '{params.dorado} basecaller '
        '{params.ext_param} '
        '-o {input.raw_bam_dir} '
        '--reference {params.reference} '
        '--mm2-opts "{params.mm2_opts}" '
        '{params.model} {input.pod5} '
        '2>>{log.err} && touch {output.tag}'

rule bam_merge:
    input:
        raw_bam_dir="data/tmp.{sample}.bam" ,
        tag="data/tmp.{sample}.bam.tag"
    output:
        bam="{project}/data/{sample}/bam/pass.bam"
    threads: config['threads']['samtools']
    conda:
        '../envs/samtools.yaml'
    priority: 10
    log:
        log="logs/bam_merge/{sample}_{project}.log",
    benchmark:
        "benchmarks/bam_merge/{sample}_{project}.txt"
    shell:
        'samtools merge '
        '-o {output.bam} '
        '--threads {threads} '
        '{input.raw_bam_dir}/*.bam'