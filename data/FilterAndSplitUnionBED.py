# Packages needed
import pandas as pd
#import gffutils
from pybedtools import BedTool
#import matplotlib.pyplot as plt
#import numpy as np
import sys
import os
#import csv

# Variables
#results_folder = "/mnt/raid/users/pilgram/methylome/analysis/processed_files/unionbedgraphs/" # Folder containing bedGraph file
#bed_file = "/mnt/raid/users/pilgram/methylome/references/targets/all_Gene_and_CpG_ID_annotations_covered_targets_Twist_Methylome_V1_TE-96341190_hg38.bed" # Folder containing bed file
#bedGraph = "bedgraphs_concat_min5_duplications.bedGraph"  # concatenated bedGraph file which contains all samples needed
#bedGraph_loc = results_folder + bedGraph 
#bedfile = BedTool(bed_file) # Load bed file as bedTool object
#sampleA = "3005027" # Sample which is going to be compared against a reference group (needed for filtering strategie 3)
#refs = ["3029134","3027101","3028523","3029136","3004931","3025912","3004928"] # Reference samples (needed for filtering strategie 3)


bedGraph = sys.argv[1]
bed_file = sys.argv[2]

outP = sys.argv[3]
bedfile = BedTool(bed_file) # Load bed file as bedTool object


bedgraph_union = pd.read_csv(bedGraph, sep="\t") # load bedGraph file into dataframe


bedgraph_strategie_1 = bedgraph_union.dropna(axis=0) # Only keep sites which are found in all samples


# Dataframe is separated to get one dataframe for each sample
bedgraphs_strategie_1 = {}
for i in bedgraph_strategie_1.columns[3:]:
    id = i
    bedgraphs_strategie_1[id] = bedgraph_strategie_1[["chrom","start","end",id]]






#os.mkdir("results_duplications_min5_strat1") # Outdir
# Map sites to enriched regions and calculate mean and stdev
for id, bedgraph in bedgraphs_strategie_1.items():
    os.makedirs(f"{outP}/{id}", exist_ok=True)
    bedtool_strategie_1 = BedTool.from_dataframe(bedgraph).sort()
    intersect_bedGraph_bed_filtered = bedfile.map(bedtool_strategie_1, nonamecheck=True, c=4, o=["count","mean","stdev"], null="N/A") # Number of sites and mean methylation per region
    intersect_bedGraph_bed_filtered.saveas(f"{outP}/{id}/{id}.bed") # Save files in outdir