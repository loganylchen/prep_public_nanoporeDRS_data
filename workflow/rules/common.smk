import glob
import pandas as pd
import sys
from snakemake.utils import validate
from snakemake.logging import logger


samples = pd.read_csv(config['samples'], sep="\t", dtype={"SampleName": str}).set_index("SampleName", drop=False).sort_index().T.to_dict()

def get_data_url(wildcards):
    return samples[wildcards.sample]["url"]




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


