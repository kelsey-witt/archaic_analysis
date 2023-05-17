"""
This script was written by Kelsey Witt Dillon in November 2020. The goal of this script is to 
examine archaic allele sharing between Europeans, East Asians, and South Asians for the purpose 
of producing a Venn Diagram-like figure. 

The input file is "archaic_rareAFR_snps.csv.gz" (infile, line 35), a csv file containing putative 
archaic alleles and a binary score for whether each Non-African individuals in the 1000 Genomes 
populations has the archaic allele. The file is generated using the script ID_rare_AFR_archaic_snps.py.

Usage: python3 archaic_site_pattern_venn.py [archaic_set], where archaic_set specifies which archaic 
individual(s) has/have the allele for it to be considered as "archaic". Options are "na_only" (found 
only in Altai Neanderthal), "nc_only" (found only in Chagyrskaya Neanderthal), "nv_only" (found only 
in Vindija Neanderthal), "d_only" (found only in Denisovan), "n_only" (found in any of the Neanderthals 
but not the Denisovan), "n_shared" (found in all three Neanderthals but not the Denisovan), "n_all" 
(found in any of the Neanderthals, regardless of presence/absence in the Denisovan), "d_all" (found 
in the Denisovan, regardless of presence/absence in any Neanderthals), "arch_shared" (shared across 
all four archaic humans), and "arch_all" (any archaic allele present in any individual).

This script generates a csv file, currently named "nonAMR_venn_counts_[archaic_set].txt" (line 36). 
Each row represents a pattern of archaic allele sharing, represented by a four digit binary value. 
The numbers represent presence or absence of the archaic allele in East Asians, Europeans, and 
South Asians (in order). A 0 means the allele is not found at least 1% frequency in a population from 
that superpopulation (according to the 1000 Genomes nomenclature), while a 1 means it is present with 
at least 1% frequency. For example, "111" means that the allele is found in all regions, while "011" 
means that the allele is present in Europeans and South Asians, but not East Asians. The second 
column of the file is the count of SNPs that show a given pattern. These numbers can be put into a 
program like venn_euler in R to make a Venn Diagram.
"""

import gzip
import sys

archaic_set = sys.argv[1] #options are na_only, nc_only, nv_only, n_shared, n_only, n_all, d_all, d_only, arch_shared, and arch_all

infile = "./archaic_rareAFR_snps.csv.gz"
outfile = "nonAMR_venn_counts_" + archaic_set + ".txt"

superPops = {"EUR": ["CEU", "TSI", "FIN", "GBR", "IBS"], \
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"], "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"]}

spops = ["EAS", "EUR", "SAS"]

superPopTracker = {}
superPopPatterns = {}

for superPop in superPops:
    superPopTracker[superPop] = []

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

g=open(outfile, 'w')

with gzip.open(infile, 'rt') as f:
    for line in f:
        if "Chromosome" in line:
            spline = line[:-1].split(sep=",")
            for colnum in range(9,len(spline)):
                indPop = spline[colnum][0:3]
                for superPop in superPops:
                    if indPop in superPops[superPop]:
                        superPopTracker[superPop].append(colnum)
        else:
            superPopHas = ""
            spline = line[:-1].split(sep=",")
            archAllele, anea, cnea, vnea, den = spline[4:9]
            archGenos = [anea, cnea, vnea, den]
            if correctArchSet: #checks to see if archaic allele is in correct individuals
                for superPop in spops:
                    inSuperPop = False
                    for col in superPopTracker[superPop]:
                        if archAllele in spline[col]:
                            inSuperPop = True
                    superPopHas += ("1" if inSuperPop else "0")
                if superPopHas in superPopPatterns:
                    superPopPatterns[superPopHas] += 1
                else:
                    superPopPatterns[superPopHas] = 1
            
spops_line = "Superpops: " + ", ".join(spops) + "\n"
g.write(spops_line)
for i in superPopPatterns:
    outLine = i + ": " + str(superPopPatterns[i]) + "\n"
    g.write(outLine)
