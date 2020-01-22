#! /usr/bin/env python

import os
from pypiper import NGSTk
import textwrap
import requests
import pandas as pd


global output_dir

SRA_URL = "https://www.ncbi.nlm.nih.gov/sra?term={url}"


EXPERIMENTS = {
    "two_level_experiment": {
        "description": "HEK293T, HeLa S3 and NIH/3T3 cell mixture (384 x 384 sci-RNA-seq)",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2599699",
        "GSE": "GSE98561",
        "GSM": "GSM2599699",
        "SRA": "SRX2784959"},
    "mid_two_level_experiment": {
        "description": "HEK293T and NIH/3T3 cell mixture (16 x 84 sci-RNA-seq)",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2599703",
        "GSE": "GSE98561",
        "GSM": "GSM2599703",
        "SRA": "SRX2784963"},
    "three_level_experiment": {
        "description": "HEK293T and NIH/3T3 cell mixture (sciRNA-seq with three levels indexing)",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2599705",
        "GSE": "GSE98561",
        "GSM": "GSM2599705",
        "SRA": "SRX2784965"},
    "sci-plex": {
        # https://science.sciencemag.org/content/suppl/2019/12/04/science.aax6234.DC1
        "description": "HEK293T and NIH/3T3 cell mixture (sci-plex)",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4150376",
        "GSE": "GSE139944",
        "GSM": "GSM4150376",
        "SRA": "SRX7101186"},
    "splitseq_3000": {
        # http://www.sciencemag.org/content/360/6385/176/suppl/DC1
        "description": "3000 cell species mixing",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3017262",
        "GSE": "GSE110823",
        "GSM": "GSM3017262",
        "SRA": "SRX3722699",
        "SRR": "SRR6750056",
    },
    "splitseq_300": {
        "description": "300 cell species mixing",
        "URL": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3017263",
        "GSE": "GSE110823",
        "GSM": "GSM3017263",
        "SRA": "SRX3722700",
        "SRR": "SRR6750057",
    },
}


def main():
    # Get the metadata and download the files
    for experiment in EXPERIMENTS.keys():
        print(experiment)
        output_dir = f"data/external/{EXPERIMENTS[experiment]['GSE']}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        r = requests.get(SRA_URL.format(url=EXPERIMENTS[experiment]["SRA"]))

        if r.ok:
            sra_table = pd.read_html(r.content.decode())[0]
        else:
            print("Failed getting run table")
            continue
        sra_table.to_csv(
            f"{output_dir}/{EXPERIMENTS[experiment]['GSM']}.annotation.csv",
            index=False)

        for sra_run in sra_table['Run']:
            cmd = f"sam-dump {sra_run} | samtools view -bS - > {output_dir}/{sra_run}.bam"
            submit_job(cmd, sra_run, output_dir, task_name="download_sra")

    # Parse the two line per paired read into single read with tags
    for experiment in EXPERIMENTS.keys():
        output_dir = f"data/external/{EXPERIMENTS[experiment]['GSE']}"
        sra_table = pd.read_csv(f"{output_dir}/{EXPERIMENTS[experiment]['GSM']}.annotation.csv")

        splitseq = experiment.startswith("splitseq")
        runnable = 'sci-rna_parse.py' if not splitseq else 'splitseq_parse.py'
        output_bam = f"{output_dir}/{sra_run}.annotated.bam" if not splitseq else f"{output_dir}/{experiment}.annotated.bam"
        for sra_run in sra_table['Run']:
            cmd = f"python3 -u /home/arendeiro/sci-rna/src/method_comparisons/{runnable} {output_dir}/{sra_run}.bam {output_bam}"
            submit_job(cmd, sra_run, output_dir, task_name="convert_scirna_to_cemm")

    # # Merge the BAM files for joint processing (not actually required but sometimes convenient)
    # for experiment in EXPERIMENTS.keys():
    #     output_dir = f"data/external/{EXPERIMENTS[experiment]['GSE']}"
    #     gsm = EXPERIMENTS[experiment]['GSM']
    #     sra_table = pd.read_csv(f"{output_dir}/{gsm}.annotation.csv")

    #     if not os.path.exists(f"data/external/{gsm}"):
    #         os.makedirs(f"data/external/{gsm}")

    #     cmd = f"sambamba merge -t 12 {output_dir}/{gsm}.annotated.bam {' '.join([f'{output_dir}/{x}.annotated.bam' for x in sra_table['Run']])}"
    #     submit_job(cmd, gsm, output_dir, task_name="merge_runs_to_experiment", queue="longq", cpus=8, mem=24000)


def submit_job(
        cmd, sra_run, output_dir,
        task_name="download_sra",
        cpus=2, mem=8000, queue="shortq"):
    tk = NGSTk()

    job_name = f"{task_name}_{sra_run}"
    job_file = f"{output_dir}/{job_name}.sh"
    log_file = f"{output_dir}/{job_name}.log"

    cmd = " " + tk.slurm_header(
        job_name, log_file, queue=queue, cpus_per_task=cpus, mem_per_cpu=mem,
    ) + cmd + ";" + tk.slurm_footer() + "\n"

    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))
    tk.slurm_submit_job(job_file)
