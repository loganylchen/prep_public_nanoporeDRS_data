import pandas as pd
import sys


samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"SampleName": str}, comment="#")
    .set_index("SampleName", drop=False)
    .sort_index()
    .T.to_dict()
)


def get_data_url(wildcards):
    url = samples[wildcards.sample]["url"]
    if url.startswith("https://") | url.startswith("http://"):
        return url
    else:
        return f"https://{url}"


def get_uncompress_command(wildcards):
    if samples[wildcards.sample]["type"] == "tar.gz":
        return "tar --ignore-failed-read --no-same-permissions  -zxvf "
    elif samples[wildcards.sample]["type"] == "tar":
        return "tar --ignore-failed-read --no-same-permissions  -xvf "
    else:
        print(f"Unknown data type for sample {wildcards.sample}")
        sys.exit(1)


def get_output_list_for_one_sample(sample):
    project = samples[sample]['project']
    return [
        f"{project}/data/{sample}/blow5/nanopore.blow5",
        f"{project}/data/{sample}/fastq/pass.fq.gz",
        f"{project}/data/{sample}/bam/pass.bam",
        f"{project}/data/{sample}/pod5/nanopore.pod5",
    ]


def get_final_output():
    final_output = []
    for sample in samples.keys():
        final_output += get_output_list_for_one_sample(sample)
    return final_output
