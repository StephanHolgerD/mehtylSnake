import pandas as pd
import os
from glob import glob

annotation_file = '/mnt/raid/users/pilgram/methylome/analysis/Homo_sapiens.GRCh38.108.chr.chr.gff3'
annotation_file_db = '/mnt/raid/users/pilgram/methylome/analysis/Homo_sapiens.GRCh38.108.chr.chr.db'
minAlignments = [5,10,20,50]

bed = 'data/all_Gene_and_CpG_ID_annotations_covered_targets_Twist_Methylome_V1_TE-96341190_hg38.bed'
df=pd.read_csv("samples.tsv",sep="\t")

SAMPLES=list(df["ID"])
SAMPLES_PATH=list(df["Bedgraph"])


def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

for a,b in zip(SAMPLES,SAMPLES_PATH):
    a=str(a)
    os.makedirs("../01_bedgraphs/"+a, exist_ok=True)
    check_symlink(b,f'../01_bedgraphs/{a}/{a}.bedgraph')



rule all:
    input:
        expand('../02_FilteredBedgraphs/{sample}/{sample}_min{filter}.bedgraph',sample=SAMPLES,filter=minAlignments),
        expand('../02_FilteredBedgraphs/{sample}/min{filter}/summary_min{filter}.csv',sample=SAMPLES,filter=minAlignments),
        #'../03_ConcatFilteredBedgraphs/ConcatFilteredBedgraphs.bed',
        #expand('../04_FilterSplitUnionBed/{sample}/{sample}.bed',sample=SAMPLES)


rule FilterBedgraphs:
    input:
        bedgraph = "../01_bedgraphs/{sample}/{sample}.bedgraph"
    threads: 1
    output:
        bedgraph_filtered = "../02_FilteredBedgraphs/{sample}/{sample}_min{filter}.bedgraph"
    shell:
        "awk '$5+$6>={wildcards.filter}' {input.bedgraph} | cut -f1-4 > {output.bedgraph_filtered}"



rule MethSeqOverview:
    input:
        bedgraph = "../01_bedgraphs/{sample}/{sample}.bedgraph"
    threads: 1
    output:
        summary = "../02_FilteredBedgraphs/{sample}/min{filter}/summary_min{filter}.csv"
    conda:
        'envs/OverviewRule.yaml'
    threads: 1
    shell:
        'python data/methylome/methylome_summary.py \
            {input.bedgraph} \
            {annotation_file} \
            {annotation_file_db} \
            {bed} \
            -min_align {wildcards.filter} \
            -outdir ../02_FilteredBedgraphs/{wildcards.sample}/min{wildcards.filter}/ \
            -outfile {output.summary}'



rule ConcatFilteredBedgraphs:
    input:
        bedgraph_filtered = expand("../02_FilteredBedgraphs/{sample}/{sample}.bedgraph",sample=SAMPLES)

    output:
        ConcatBed = '../03_ConcatFilteredBedgraphs/ConcatFilteredBedgraphs.bed'

    threads: 1

    conda:
        'envs/bedtools.yaml'

    shell:
        'bedtools unionbedg -i {input.bedgraph_filtered} -header -names {SAMPLES} -filler N/A > {output.ConcatBed}'



rule FilterAndSplitUnioBed:
    input:
        ConcatBed = '../03_ConcatFilteredBedgraphs/ConcatFilteredBedgraphs.bed'

    output:
        SplitBed = expand('../04_FilterSplitUnionBed/{sample}/{sample}.bed',sample=SAMPLES)

    conda:
        'envs/pybedtools.yaml'

    threads: 1

    shell:
        'python data/FilterAndSplitUnionBED.py {input.ConcatBed} {bed} ../04_FilterSplitUnionBed/'
