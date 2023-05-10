"""
This script was written by Kelsey Witt Dillon in June 2021. The goal of this script is to 
count the number of archaic alleles found in each non-African individual in the 1000 Genomes
populations, and to also calculate allele frequencies for both archaic and modern (not shared
with archaic genomes) alleles with low frequency (<1%) in Africa, and classify them as rare 
(above 1% frequency but below 20% frequency) or common (20% frequency or greater). For an allele 
to be classified as archaic, it had to be found in Africans at a frequency less than 1% but at 
a frequency of 1% or more in at least one non-African population and shared with the archaic 
humans of interest (as established by archaic_set).

The input files are the csv files for each chromosome containing alleles that are rare in Africa 
(which come from ID_rare_AFR_genotypes.py, "afrInfile", line 160), and the four archaic human VCF
files ("aNeaInfile", "vNeaInfile","cNeaInfile", "denInfile", lines 161-164.)

Usage: python3 rareAFR_allele_counter.py [archaic_set], where archaic_set specifies which archaic individual(s) 
has/have the allele for it to be considered as "archaic". Options are "na_only" (found only in Altai Neanderthal),
"nc_only" (found only in Chagyrskaya Neanderthal), "nv_only" (found only in Vindija Neanderthal), "d_only"
(found only in Denisovan), "n_only" (found in any of the Neanderthals but not the Denisovan), "n_shared"
(found in all three Neanderthals but not the Denisovan), "n_all" (found in any of the Neanderthals, regardless
of presence/absence in the Denisovan), "d_all" (found in the Denisovan, regardless of presence/absence in any
Neanderthals), "arch_shared" (shared across all four archaic humans), and "arch_all" (any archaic allele
present in any individual).

This script generates two output files. The first, "[archaic_set]_arch_allele_counts.csv" (popOutfile, line 45) 
is a csv file that summarizes population-level information. Each row is a population, and the columns are, 
in order, the population name, the number of non-archaic alleles that are rare in Africa but 20%+ frequency 
in that population ("modern common"), the number of non-archaic alleles that are rare in Africa but between 
1 and 20% frequency in that population ("modern rare"), the number of archaic alleles that are rare in Africa 
but 20%+ frequency in that population ("archaic common"), and the number of archaic alleles that are rare in 
Africa but between 1 and 20% frequency in that population ("archaic rare"). 

The second output file, "[archaic_set]_ind_arch_score_singles.csv" (indOutfile, line 46) is a csv file that 
summarizes individual-level information. Each row represents an individual, and the columns are, in order, the 
individual ID, the population the individual is from, and the number of archaic alleles that individual has. 
"Archaic" alleles are determined by the criteria above, and are counted across all sites in the genome (ie 
a site that is homozygous for an archaic allele would increase the total by 2).
"""

import gzip
import csv
import sys

archaic_set = sys.argv[1] #possible options: arch_all, arch_shared, n_all, n_only, na_only, nc_only, nv_only, n_shared, d_all, d_only

popOutfile = archaic_set + "_arch_allele_counts.csv"
indOutfile = archaic_set + "_ind_arch_score_singles.csv"

archLineList = ["aLine", "cLine", "vLine", "dLine"]
archLines = {"aLine":"", "cLine":"", "vLine":"", "dLine":""}
archPos = [0,0,0,0]
archLims = {"2": 243177959, "3":197874834, "4":191043233, "5":180715571, "6": 170919480, "7": 159128533, "8": 146300828, "9": 141102551, \
    "10": 135512040, "11": 134946370, "12": 133841510, "13": 115108863, "14": 107289205, "15": 102506248, "16": 90174509, "17": 81160199,  \
        "18": 78017036, "20": 62917942, "21": 48101478, "22": 51234514}

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))
chr_list.append("X")

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Pop_Tracker = {}
Ind_Tracker = {}
Modern_Allele_Counter = {}
Archaic_Allele_Counter = {}
Ind_Allele_Counter = {}

for pop in populations:
    Pop_Tracker[pop] = []
    Modern_Allele_Counter[pop] = {}
    Modern_Allele_Counter[pop]["common"] = 0
    Modern_Allele_Counter[pop]["rare"] = 0
    Archaic_Allele_Counter[pop] = {}
    Archaic_Allele_Counter[pop]["common"] = 0
    Archaic_Allele_Counter[pop]["rare"] = 0

def skip_archaic_header(aFile,line):
    while '#' in line:
        line = aFile.readline()
    return line

def scan_archaic_lines(aFile,archaicPos,afrPos):
    if chromosome in archLims:
        while int(archaicPos)<int(afrPos) and int(archaicPos)<archLims[chromosome]:
            line = aFile.readline()
            spline = line.split()
            archaicPos = spline[1]
    else:
        while int(archaicPos)<int(afrPos):
            line = aFile.readline()
            spline = line.split()
            archaicPos = spline[1]
    return line

def calc_total_X_allele(allele_total, genotype):
    allele_count_x = sum(genotype.count(x) for x in ("0", "1"))
    allele_total += allele_count_x
    return allele_total

def calc_pop_freq(population, nonafr_allele):
    pop_allele_count = 0
    allele_total = len(Pop_Tracker[population]) * 2
    if chr == "X":
        allele_total = 0
    for ind in Pop_Tracker[population]:
        pop_allele_count += spline[ind].count(nonafr_allele)
        if chr == "X": #count number of X alleles as we go to account for recombinant and non-recombinant
            allele_total = calc_total_X_allele(allele_total, spline[ind])
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return rounded_freq

def allele_counter(freq_list, counter_dict):
    for x in range(0,len(freq_list)):
        if float(freq_list[x]) > 0.20:
            counter_dict[populations[x]]["common"] += 1
        elif float(freq_list[x]) > 0.01:
            counter_dict[populations[x]]["rare"] += 1

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
      
for chromosome in chr_list:
    afrInfile = "rare_AFR_genotypes_derived_chr" + chromosome + ".csv.gz"
    aNeaInfile = "/users/kwittdil/data/data/archaic_genomes/altai_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    cNeaInfile = "/users/kwittdil/data/data/archaic_genomes/chagyrskaya_VCF/ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr" + chromosome + ".noRB.vcf.gz"
    vNeaInfile = "/users/kwittdil/data/data/archaic_genomes/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    denInfile = "/users/kwittdil/data/data/archaic_genomes/denisova_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    af = gzip.open(aNeaInfile, 'rt')
    cf = gzip.open(cNeaInfile, 'rt')
    vf = gzip.open(vNeaInfile, 'rt')
    df = gzip.open(denInfile, 'rt')
    archFiles = [af, cf, vf, df]
    for i in range(0, len(archFiles)):
        archaicLine = archLineList[i]
        archLines[archaicLine] = archFiles[i].readline()
        afile = archFiles[i]
        aLine = archLines[archaicLine]
        archLines[archaicLine] = skip_archaic_header(afile,aLine)
        archSpline = archLines[archaicLine].split()
        archPosition = archSpline[1]
        archPos[i] = archPosition
    with gzip.open(afrInfile, 'rt') as f:
        for line in f:
            if "Chromosome" in line and chromosome == "1":
                spline = line[:-1].split(sep=",")
                for colnum in range(5,len(spline)):
                    indID = spline[colnum]
                    Ind_Tracker[colnum] = indID
                    Ind_Allele_Counter[colnum] = 0
                    indPop = indID[0:3]
                    Pop_Tracker[indPop].append(colnum)
            elif "Chromosome" not in line:
                freq_by_pop = []
                arch_genotypes = ["-","-","-","-"]
                spline = line[:-1].split(sep=",")
                position, ref_allele, alt_allele, rare_allele = spline[1:5]
                if chromosome not in archLims or (chromosome in archLims and int(position) <= archLims[chromosome]):
                    for j in range(0, len(archFiles)):
                        archaicPosition = archPos[j]
                        archaicLine = archLineList[j]
                        currArchLine = archLines[archaicLine]
                        afile = archFiles[j]
                        if int(archaicPosition) < int(position):
                            archLines[archaicLine] = currArchLine = scan_archaic_lines(afile,archaicPosition,position)
                        aSpline = currArchLine.split()
                        archaicPosition,_,arch_ref,arch_alt,genotype_qual = aSpline[1:6]
                        archPos[j] = archaicPosition
                        if int(archaicPosition) == int(position):
                            if arch_alt in [".",alt_allele]:
                                if float(genotype_qual) >= 40:
                                    arch_genotypes[j] = aSpline[9][0:3]
                    for pop in populations:
                        nonafr_freq = calc_pop_freq(pop, rare_allele)
                        freq_by_pop.append(nonafr_freq)
                        for ind in Pop_Tracker[pop]:
                            if correctArchSet(rare_allele,arch_genotypes): 
                                Ind_Allele_Counter[ind] += spline[ind].count(rare_allele)
                    if correctArchSet(rare_allele,arch_genotypes): 
                        allele_counter(freq_by_pop, Archaic_Allele_Counter)
                    else:
                        allele_counter(freq_by_pop, Modern_Allele_Counter)

with open(indOutfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["ind", "pop", "positions"]
    w.writerow(header)
    for pop in populations:
        for ind in Pop_Tracker[pop]:
            indID = Ind_Tracker[ind]
            indScore = str(Ind_Allele_Counter[ind])
            score_line = [indID,pop,indScore]
            w.writerow(score_line)

with open(popOutfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["pop", "mod_common", "mod_rare", "archaic_common", "archaic_rare"]
    w.writerow(header)
    for i in range(0,len(populations)):
        allele_line = [populations[i], str(Modern_Allele_Counter[populations[i]]["common"]), str(Modern_Allele_Counter[populations[i]]["rare"]),  str(Archaic_Allele_Counter[populations[i]]["common"]), str(Archaic_Allele_Counter[populations[i]]["rare"])]
        w.writerow(allele_line)
