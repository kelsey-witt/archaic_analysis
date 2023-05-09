"""
This script was written by Kelsey Witt Dillon in May 2021. The goal of this script is to 
identify the number of archaic alleles in a subset of individuals from each population in
the 1000 Genomes Project. Sample sizes range from 1-150, and the script samples 100 replicates
for each sample size. 

The input file is "archaic_rareAFR_snps.csv.gz" (infile, line 27), a csv file containing putative archaic alleles
and a binary score for whether each Non-African individuals in the 1000 Genomes populations has
the archaic allele. The file is generated using the script ID_rare_AFR_archaic_snps.py. The list of sample_sizes
is set to N=1, 5, 10, 25, 50, 75, 100, 125, and 150, but can be adjusted in line 36.

Usage: python3 archaic_mp_cvg_ss.py [archaic_set], where archaic_set specifies which archaic individual(s) 
has/have the allele for it to be considered as "archaic". Options are "na_only" (found only in Altai Neanderthal),
"nc_only" (found only in Chagyrskaya Neanderthal), "nv_only" (found only in Vindija Neanderthal), "d_only"
(found only in Denisovan), "n_only" (found in any of the Neanderthals but not the Denisovan), "n_shared"
(found in all three Neanderthals but not the Denisovan), "n_all" (found in any of the Neanderthals, regardless
of presence/absence in the Denisovan), "d_all" (found in the Denisovan, regardless of presence/absence in any
Neanderthals), "arch_shared" (shared across all four archaic humans), and "arch_all" (any archaic allele
present in any individual).

This script generates a csv file for each sample size, labeled "[archaic_set]_mp_cvg_[N]ind_100x.csv", where N
is the sample size, and archaic_set is the variable specified above. Each row is a replicate and each column is
a population, and the number of unique archaic sites found across the sample size of that number of individuals
from a given population is listed in the table.
"""

import csv
import random
import sys
import gzip

archaic_set = sys.argv[1] #possible options: arch_all, arch_shared, n_all, n_only, na_only, nc_only, nv_only, n_shared, d_all, d_only

infile = "./archaic_rareAFR_snps.csv.gz"
sample_sizes = [1,5,10,25,50,75,100,125,150,170]

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Pop_Tracker = {}
Pop_Cvg = {}

def correctArchSet(allele, genotypes): 
    fitsCriteria = False
    anea, cnea, vnea, den = genotypes[0:4]
    if archaic_set == "arch_all":
        fitsCriteria = any(allele in geno for geno in genotypes)
    elif archaic_set in ["d_only", "na_only", "nc_only", "nv_only", "n_only", "n_shared", "arch_shared"]:
        if all(geno != "-" for geno in genotypes):
            if archaic_set == "d_only":
                if allele in den and all(allele not in geno for geno in [anea, cnea, vnea]):
                    fitsCriteria = True
            elif archaic_set == "na_only":
                if allele in anea and all(allele not in geno for geno in [den, cnea, vnea]):
                    fitsCriteria = True
            elif archaic_set == "nc_only":
                if allele in cnea and all(allele not in geno for geno in [den, anea, vnea]):
                    fitsCriteria = True
            elif archaic_set == "nv_only":
                if allele in vnea and all(allele not in geno for geno in [den, anea, cnea]):
                    fitsCriteria = True
            elif archaic_set == "n_only":
                if allele not in den and any(allele in geno for geno in [anea, cnea, vnea]):
                    fitsCriteria = True        
            elif archaic_set == "n_shared":
                if allele not in den and all(allele in geno for geno in [anea, cnea, vnea]):
                    fitsCriteria = True    
            elif archaic_set == "arch_shared":
                if all(allele in geno for geno in genotypes):
                    fitsCriteria = True
    elif archaic_set == "d_all":
        if allele in den:
            fitsCriteria = True
    elif archaic_set == "n_all":
        fitsCriteria = any(allele in geno for geno in [anea, cnea, vnea])
    else:
        print("You have selected a non-standard archaic set")
    if fitsCriteria:
        return True
    else:
        return False

for ss in sample_sizes:
    outfile = archaic_set + "_mp_cvg_" + str(ss) + "ind_100x.csv"
    with open(outfile, 'w', newline = '') as csvfile:
        w = csv.writer(csvfile, delimiter=",")
        header = ["run"] + populations
        w.writerow(header)
        for i in range(1,10):
            for pop in populations:
                Pop_Cvg[pop] = 0
                Pop_Tracker[pop] = []
            with gzip.open(infile, 'rt') as f:
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
                        archAllele = spline[4]
                        archGenos = spline[5:9]
                        if correctArchSet(archAllele,archGenos):
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
