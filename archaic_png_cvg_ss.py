"""
This script was written by Kelsey Witt Dillon in November 2021. The goal of this script is to 
identify the number of archaic alleles in a subset of individuals from each population in
the 1000 Genomes Project and Papuans from the Simons Genome Diversity Project. Sample sizes 
range from 1-15, and the script samples 100 replicates for each sample size. 

The input file is "[archaic_set]_snps_per_ind_png.csv" (infile, line 35), a csv file containing 
putative archaic alleles and a binary score for whether each Non-African individual in the 1000 
Genomes populations and each Papuan from the Simons Genome Diversity Project has the archaic allele. 
The file is generated using the script allelecvg_make_table_v2_png.py. The list of sample_sizes
is set to N=1, 5, 10, and 15, (ie the total number of Papuan genomes) but can be adjusted in line 36.

Usage: python3 archaic_mp_cvg_ss.py [archaic_set], where archaic_set specifies which archaic individual(s) 
has/have the allele for it to be considered as "archaic". Options are "d_only" (found only in Denisovan), 
"n_only" (found in any of the Neanderthals but not the Denisovan), and "arch_all" (any archaic allele
present in any individual).

This script generates a csv file for each sample size, labeled "[archaic_set]_png_cvg_[N]ind_100x.csv", where N
is the sample size, and archaic_set is the variable specified above. Each row is a replicate and each column is
a population, and the number of unique archaic sites found across the sample size of that number of individuals
from a given population is listed in the table.
"""

import csv
import random
import sys
import gzip

archaic_set = sys.argv[1] #possible options: arch_all, n_only, d_only

infile = "./" + archaic_set + "_snps_per_ind_png.csv"
sample_sizes = [1,5,10,15]

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU","PNG"]
Pop_Tracker = {}
Pop_Cvg = {}

for ss in sample_sizes:
    outfile = archaic_set + "_png_cvg_" + str(ss) + "ind_100x.csv"
    with open(outfile, 'w', newline = '') as csvfile:
        w = csv.writer(csvfile, delimiter=",")
        header = ["run"] + populations
        w.writerow(header)
        for i in range(1,100):
            for pop in populations:
                Pop_Cvg[pop] = 0
                Pop_Tracker[pop] = []
            with open(infile) as f:
                for line in f:
                    if 'Chromosome' in line:
                        spline = line[:-1].split(sep=",")
                        for colnum in range(9,len(spline)):
                            indPop = spline[colnum][0:3]
                            Pop_Tracker[indPop].append(colnum)
                        for pop in Pop_Tracker.keys():
                            ind_downsample = random.choices(Pop_Tracker[pop], k=ss)
                            ind_downsample.sort()
                            Pop_Tracker[pop] = ind_downsample
                    else:
                        spline = line[:-1].split(sep=",")
                        for pop in populations:
                            hasArchAllele = False
                            for ind in Pop_Tracker[pop]:
                                if "1" in spline[ind]:
                                    hasArchAllele = True
                            if hasArchAllele:
                                Pop_Cvg[pop] += 1
            outline = [str(i)]
            for pop in populations:
                outline.append(str(Pop_Cvg[pop]))
            w.writerow(outline)
