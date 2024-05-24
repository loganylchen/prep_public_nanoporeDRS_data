#!/usr/bin/env bash
set -x
set -e



lftp -c "pget -n ${snakemake[threads]} ${snakemake[params]} -o {output.compressed_data}" >{log}
touch {output.tag}