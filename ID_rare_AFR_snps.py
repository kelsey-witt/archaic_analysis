"""
This script was written by Kelsey Witt Dillon in April 2021. The goal of this script is to 
make a file containing all biallelic sites from the 1000 Genomes dataset where one allele has a 
frequency less than 1% in African populations and a minimum frequency of 1% in at least one non-African 
population. It records whether or not each non-African individual has the allele that is rare in 
Africa at that position.

The input files are the .panel file from the 1000 Genomes dataset ("pop_file", line 27), and the VCF
files of each chromosome from the 1000 Genomes Phase III Dataset ("modern_file", lines 110 and 112).

Usage: python3 ID_rare_AFR_snps.py

This script generates a csv file for each chromosome, currently named "rare_AFR_snps_chrN.csv", where
N is 1-22 or X ("outfile", line 114). Each row represents a site in the genome, and the first five columns
are the chromosome number, position, the reference allele, the alternate allele, and the allele that is rare
in Africa (0 if reference, 1 if alternate). The remaining columns are the presence or absence of the rare allele
in that individual's genotype, organized by population. A "0" means that the individual lacks the non-African 
allele, while a "1" indicates that the non-African allele is present in that individual. The header lists the 
individual's ID as their population and the 0-indexed column number they appear in the 1000 Genomes vcf, 
separated by an underscore (ie "GBR_9").
"""

import gzip
import re
import csv

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
afrPops = ["YRI", "LWK", "GWD", "MSL", "ESN"]
Pop_Tracker = {}

for pop in populations:
    Pop_Tracker[pop] = []

chr_list = ["X"]
for c in range(1,23):
    chr_list.append(str(c))

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] == "AFR" and line_col[1] in afrPops: 
            Pop_Tracker["AFR"].append(col_counter)
        elif line_col[2] != "AFR":
            pop = line_col[1]
            Pop_Tracker[pop].append(col_counter)
        col_counter += 1

def decode_split_line(line):
    #line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line[:-1].split()
    return line

def snp_is_biallelic(line):
    ref_allele = line[3]
    alt_allele = line[4]
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return True
    else:
        return False

def calc_total_X_allele(allele_total, genotype):
    allele_count_x = sum(genotype.count(x) for x in ("0", "1"))
    allele_total += allele_count_x
    return allele_total

def parse_mod_line(line):
    position, rsVal, ref_allele, mod_alt_allele = line[1:5]
    mod_info = line[7]
    return [position, ref_allele, mod_alt_allele, mod_info, rsVal]

def calc_afr_allele_freq(line):
    afr_allele_total = len(Pop_Tracker["AFR"]) * 2
    if chr == "X":
        afr_allele_total = 0
    afr_allele_1_sum = 0
    for afr_ind in Pop_Tracker["AFR"]:
        afr_allele_1_sum += line[afr_ind].count("1")
        if chr == "X": #count number of X alleles as we go to account for recombinant and non-recombinant
            afr_allele_total = calc_total_X_allele(afr_allele_total, line[afr_ind])
    afr_allele_1_fq = float(afr_allele_1_sum/afr_allele_total)
    afr_allele_1_round = round(afr_allele_1_fq,3)
    return afr_allele_1_round

def identify_non_afr_alleles(line):
    non_afr_geno = "9"
    afr_alt_af = calc_afr_allele_freq(line)
    if float(afr_alt_af) < 0.01 or float(afr_alt_af) > 0.99:
        if float(afr_alt_af) < 0.01: #allele of interest is alternate
            non_afr_geno = "1"
        elif float(afr_alt_af) > 0.99: #allele of interest is reference
            non_afr_geno = "0"
    return non_afr_geno

def calc_pop_freq(population, nonafr_allele):
    pop_allele_count = 0
    allele_total = len(Pop_Tracker[population]) * 2
    for ind in Pop_Tracker[population]:
        pop_allele_count += mod_line[ind].count(nonafr_allele)
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return rounded_freq

for chromosome in chr_list:
    #modern_file = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    modern_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    if chromosome == "X":
        modern_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    with gzip.open(modern_file, 'rt') as mf:
        outfile = "rare_AFR_snps_chr" + chromosome + ".csv"
        with open(outfile, 'w', newline = '') as csvfile:
            w = csv.writer(csvfile, delimiter=",")
            header = ["Chromosome","Position","Reference","Alternate","Rare Allele"]
            for pop in populations[1:]:
                for ind in Pop_Tracker[pop]:
                    indTitle = pop + "_" + str(ind)
                    header.append(indTitle)
            w.writerow(header)
            for mod_line in mf:
                mod_line = decode_split_line(mod_line)
                #print(mod_line)
                if '#' not in mod_line:
                    mod_vars = parse_mod_line(mod_line)
                    mod_position, mod_ref, mod_alt = mod_vars[0:3]
                    freq_by_pop = []
                    if snp_is_biallelic(mod_line):
                        for pop in populations[1:]:
                                alleleCounter = 0
                                allele_total = len(Pop_Tracker[pop]) * 2
                                if chromosome == "X":
                                    allele_total=0
                                    for ind in Pop_Tracker[pop]:
                                        allele_total = calc_total_X_allele(allele_total, mod_line[ind])
                                for ind in Pop_Tracker[pop]:
                                    alleleCounter += mod_line[ind].count("1")
                        non_afr_genotype = identify_non_afr_alleles(mod_line)
                        if non_afr_genotype == "0" or non_afr_genotype == "1": #if there is an allele that's rare in Africans
                            nonRare = False
                            for pop in populations:
                                nonafr_freq = calc_pop_freq(pop, non_afr_genotype)
                                freq_by_pop.append(str(nonafr_freq))
                                if pop != "AFR" and float(nonafr_freq) > 0.01:
                                    nonRare = True
                            if float(freq_by_pop[0])<0.01 and nonRare:
                                outLine = [chromosome, mod_position, mod_ref, mod_alt, non_afr_genotype]
                                for pop in populations[1:]:
                                    nonRareInPop = False
                                    nonafr_freq = calc_pop_freq(pop, non_afr_genotype)
                                    if nonafr_freq > 0.01: #checks for frequency in individual populations
                                        nonRareInPop = True 
                                    for ind in Pop_Tracker[pop]:
                                        if non_afr_genotype in mod_line[ind] and nonRareInPop:
                                            indVal = "1"
                                        else:
                                            indVal = "0"
                                        outLine.append(indVal)
                                w.writerow(outLine)