#!/usr/bin/env python

import os
from ngs_toolkit.utils import download_gzip_file

# initialize project
root_output_dir = "/home/arendeiro/resources/genomes"


organisms = ['homo_sapiens', 'drosophila_melanogaster']
genome_provider = ['ensembl', 'ensembl']
genome_provider_releases = ['97', '97']
assemblies = ['GRCh38', 'BDGP6.22']

genome_handle = "_".join(assemblies)
genome_dir = os.path.join(root_output_dir, genome_handle)
genomes_dir = [os.path.join(genome_dir, assembl) for assembl in assemblies]

for _dir in [root_output_dir, genome_dir] + genomes_dir:
    if not os.path.exists(_dir):
        os.makedirs(_dir)

provider_url = {
    "ensembl": "ftp://ftp.ensembl.org/pub/"
}

args = zip(organisms, genome_provider, genome_provider_releases, assemblies, genomes_dir)

for organism, genome_provider, release, assembly, genome_dir in args:
    genome_file = f"{organism.capitalize()}.{assembly}.dna.toplevel.fa.gz"
    url = f"{provider_url[genome_provider]}release-{release}/fasta/{organism}"
    url += f"/dna_index/{genome_file}"
    download_gzip_file(url, os.path.join(genome_dir, genome_file))


    # Add extra chromosomes (constructs) to genome
    cmd = (
        "cat {} {} > {}"
        .format(
            os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
            output_fasta,
            os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa")))
    os.system(cmd)
    cmd = (
        "cat {} {} > {}"
        .format(
            os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"),
            output_gtf,
            os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf")))
    os.system(cmd)

    # Build STAR index (contruct + spiked with gRNAs)
    cmd = "srun --mem 80000 -p develop /cm/shared/apps/star/2.4.2a/STAR"
    cmd += " --runThreadN 8"
    cmd += " --runMode genomeGenerate"
    cmd += " --genomeDir {}".format(spiked_dir)
    cmd += " --genomeFastaFiles {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
    cmd += " --sjdbGTFfile {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
    cmd += " --sjdbOverhang 74"
    os.system(cmd)

    # Create sequence dictionaries (for piccard)
    cmd = "srun --mem 80000 -p develop java -Xmx8g -jar /cm/shared/apps/picard-tools/1.140/picard.jar"
    cmd += " CreateSequenceDictionary"
    cmd += " REFERENCE={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
    cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
    cmd += " GENOME_ASSEMBLY={}".format(genome_ref)
    cmd += " SPECIES=human"
    os.system(cmd)

    # Create reflat files
    cmd = "srun --mem 80000 -p develop java -Xmx80g -jar ~/Drop-seq_tools-1.12/jar/dropseq.jar ConvertToRefFlat"
    cmd += " SEQUENCE_DICTIONARY={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.dict"))
    cmd += " ANNOTATIONS_FILE= {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
    cmd += " OUTPUT={}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.refFlat"))
    os.system(cmd)

    # Remove vanilla genome
    os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa"))
    os.remove(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.gtf"))


    # Build STAR index for newer version (contruct + spiked with gRNAs)
    cmd = "srun -J new_STAR_reference --mem 80000 -p develop /workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR"
    cmd += " --runThreadN 8"
    cmd += " --runMode genomeGenerate"
    cmd += " --genomeDir {}".format(spiked_dir)
    cmd += " --genomeFastaFiles {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa"))
    cmd += " --sjdbGTFfile {}".format(os.path.join(spiked_dir, "Homo_sapiens.GRCh38.77.spiked.gtf"))
    cmd += " --sjdbOverhang 74"
    os.system(cmd)

    # sbatch \
    # -J new_STAR_reference \
    # -c 12 --mem 80000 -p develop \
    # --wrap "/cm/shared/apps/star/2.4.2a/STAR \
    # --runThreadN 8 \
    # --runMode genomeGenerate \
    # --genomeDir /home/arendeiro/resources/genomes/hg38_spiked_lambda/indexed_star_index_2.4.2a \
    # --genomeFastaFiles /home/arendeiro/resources/genomes/hg38_spiked_lambda/Homo_sapiens.GRCh38.dna.primary_assembly.spiked.fa \
    # --sjdbGTFfile /home/arendeiro/resources/genomes/hg38_spiked_lambda/Homo_sapiens.GRCh38.77.spiked.gtf \
    # --sjdbOverhang 74"


    # sbatch \
    # -J new_STAR_reference \
    # -c 12 --mem 80000 -p develop \
    # --wrap "~/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR \
    # --runThreadN 12 \
    # --runMode genomeGenerate \
    # --genomeDir /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e \
    # --genomeFastaFiles /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.fa \
    # --sjdbGTFfile /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf \
    # --sjdbOverhang 74"


    # sbatch -c 24 -p mediumq --time 1-08:00:00 -J new_STAR_reference --mem 80000 \
    # --wrap "~/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR \
    # --runThreadN 24 \
    # --runMode genomeGenerate \
    # --genomeDir ~/resources/genomes/hg38/indexed_STAR-2.7.0e \
    # --genomeFastaFiles /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa \
    # --sjdbGTFfile /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf \
    # --sjdbOverhang 54"
