"""
This script was written by Kelsey Witt Dillon in May 2021. The goal of this script is to 
take the csv file of archaic genotypes generates using "ID_rare_AFR_archaic_genotypes.py"
and converting it to the input files needed to run smartpca.

The input files are the .panel file for the 1000 Genomes Project (pop_file, line 31), and 
the archaic genotype file generated using ID_rare_AFR_archaic_genotypes.py (modern_file, 
line 103)

Usage: python3 cvgtable_to_smartpca.py [archaic_set], where archaic_set (line 28) specifies 
which archaic individual(s) has/have the allele for it to be considered as "archaic". 
Options are "d_only" (found only in Denisovan), "n_only" (found in any of the Neanderthals 
but not the Denisovan),  and "arch_all" (any archaic allele present in any individual).

This script generates three files, which are the required inputs for smartpca. The file 
suffix is defined by outfile_name (line 32) and there are 3 files: a .genotypes file,
a .snps file, and a .ind file. Briefly, the .ind file is a list of all individuals
in order, the .snps file details position information about the snps, and the .genotypes
file supplies genotype information for each individual in the form of the number
of reference alleles they possess.
"""

import gzip
import re
import sys
import csv

archaic_set = sys.argv[1] #possible options: d_only, n_only, arch_all

#pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
outfile_name = "mp_pca_" + archaic_set
genotype_outfile = outfile_name + ".genotypes"
snp_outfile = outfile_name + ".snps"
ind_outfile = outfile_name + ".ind"

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Pop_Tracker = {}
popTracker = {}

for pop in populations:
    Pop_Tracker[pop] = []
    popTracker[pop] = []

with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] != "AFR":
            ind,pop,_,sex = line_col[0:4]
            Pop_Tracker[pop].append([ind,sex])

with open(ind_outfile, 'w') as f:
    for pop in populations:
        for i in range(0,len(Pop_Tracker[pop])):
            indID,sex = Pop_Tracker[pop][i][0:2]
            indLine = indID + "\t" + sex + "\t" + pop + "\n"
            f.write(indLine)

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

geno = open(genotype_outfile, 'w')
snp = open(snp_outfile, 'w')

modern_file =  "archaic_rareAFR_genotypes.csv.gz"
with gzip.open(modern_file, 'rt') as mf:
    for mod_line in mf:
        if "Chromosome" in mod_line:
            line = mod_line.split(sep=",")
            for colnum in range(9,len(line)):
                indPop = line[colnum][0:3]
                popTracker[indPop].append(colnum)
        else:
            line = mod_line.split(sep=",")
            chromosome, mod_position, ref_allele, alt_allele, arch_allele = line[0:5]
            archGenos = line[5:9]
            if correctArchSet(arch_allele,archGenos):
                aboveFive = False
                for pop in populations:
                    totalChr = len(popTracker[pop])*2
                    archaicChr = 0
                    for ind in popTracker[pop]:
                        archaicChr += line[ind].count(str(arch_allele))
                    alleleFreq = float(archaicChr/totalChr)
                    if alleleFreq >= 0.05:
                        aboveFive = True
                if aboveFive:
                    if chromosome == "X":
                        chromosome = "23"
                    snpLine = ". \t" + chromosome + "\t 0.0 \t" + str(mod_position) + "\t" + ref_allele + "\t" + alt_allele + "\n"
                    snp.write(snpLine)
                    genotypeList = ""
                    for col in range(9,len(line)):
                        numRef = str(line[col].count("0"))
                        genotypeList += numRef
                    genotypeList += "\n"
                    geno.write(genotypeList)
geno.close()
snp.close()
