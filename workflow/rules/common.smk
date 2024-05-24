import pandas as pd
import sys



samples = pd.read_csv(config['samples'], sep="\t", dtype={"SampleName": str},comment='#').set_index("SampleName", drop=False).sort_index().T.to_dict()

def get_data_url(wildcards):
    return samples[wildcards.sample]["url"]


def get_compressed_data(wildcards):
    if samples[wildcards.sample]['type'] == 'tar.gz':
        return {
            'compressed_data':temp(f"data/tmp.{wildcards.sample}.nanoporeDRS.tar.gz"),
            'tag': temp(f"data/tmp.{wildcards.sample}.nanoporeDRS.download.tag")
        }
    elif samples[wildcards.sample]['type'] == 'tar':
        return {
            'compressed_data':temp(f"data/tmp.{wildcards.sample}.nanoporeDRS.tar"),
            'tag': temp(f"data/tmp.{wildcards.sample}.nanoporeDRS.download.tag")
        }
    else:
        print(f"Unknown data type for sample {wildcards.sample}")
        sys.exit(1)

def get_compressed_data_notemp(wildcards):
    if samples[wildcards.sample]['type'] == 'tar.gz':
        return {
            'compressed_data': f"data/tmp.{wildcards.sample}.nanoporeDRS.tar.gz",
            'tag': f"data/tmp.{wildcards.sample}.nanoporeDRS.download.tag"
        }
    elif samples[wildcards.sample]['type'] == 'tar':
        return {
            'compressed_data': f"data/tmp.{wildcards.sample}.nanoporeDRS.tar",
            'tag': f"data/tmp.{wildcards.sample}.nanoporeDRS.download.tag"
        }
    else:
        print(f"Unknown data type for sample {wildcards.sample}")
        sys.exit(1)

def get_uncompress_command(wildcards):
    if samples[wildcards.sample]['type'] == 'tar.gz':
        return {
            'uncompress': "tar -zxvf ",

        }
    elif samples[wildcards.sample]['type'] == 'tar':
        return {
            'uncompress': "tar -xvf ",
        }
    else:
        print(f"Unknown data type for sample {wildcards.sample}")
        sys.exit(1)


def get_output_list_for_one_sample(sample):
    return [
        # f"data/{sample}/fastq/pass.fq.gz",
        f"data/{sample}/blow5/nanopore.blow5",
    ]

def get_final_output():
    final_output = []
    for sample in samples.keys():
        final_output += get_output_list_for_one_sample(sample)
    return final_output


