import os
from snakemake.shell import shell


def extract_tar_files(directory,compressed_files=[]):

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            if file_path.endswith(".tar") or file_path.endswith(".tar.gz"):
                compressed_files.append(file_path)
            else:
                print(f"Skipping {filename}")
        elif os.path.isdir(file_path):
            compressed_files = (file_path,compressed_files)
    return compressed_files

compressed_files = []
compressed_files = extract_tar_files(snakemake.input[0],compressed_files)

for f in compressed_files:
    f_dir = os.path.dirname(f)
    if f.endswith('.tar.gz'):

        shell("tar -zxvf "+f+" -c  " + f_dir)
    else:
        shell(f"tar -xvf "+f+" -c "+f_dir)

shell("touch {snakemake.output[0]}")