import os
from snakemake.shell import shell


def extract_tar_files(directory,compressed_files=[]):

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            if filename.endswith(".tar") or filename.endswith(".tar.gz"):
                compressed_files.append(file_path)
            else:
                print(f"Skipping {filename}")
        elif os.path.isdir(file_path):
            compressed_files = (file_path,compressed_files)
    return compressed_files

compressed_files = []
compressed_files = extract_tar_files(snakemake.input[0],compressed_files)

for f in compressed_files:
    if f.endswith('.tar.gz'):
        shell("tar -zxvf {f} -c dirname({f}) ")
    else:
        shell("tar -xvf {f} -c dirname({f}) ")

shell("touch {snakemake.output[0]}")