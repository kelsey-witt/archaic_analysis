"""
This script was written by Kelsey Witt Dillon in January 2022. The goal of this script is to 
examine archaic allele sharing between Europeans, East Asians, South Asians, and Papuans
for the purpose of producing a Venn Diagram-like figure. Because the Simons Genome Diversity
Project only has 15 Papuans, the other populations are downsampled to match the same sample
size.

The input files are the .panel file and VCF files from the 1000 Genomes Phase 3 dataset (pop_file 
in line 47 and kgFile in line 243), a merged vcf containing the genotype calls for the SGDP Papuans
separated by chromosome (papuanInfile, line 248) and the four archaic human VCF files ("aNeaInfile", 
"vNeaInfile","cNeaInfile", "denInfile", lines 244-247.)

Usage: python3 archaic_venn_1kg_png.py [archaic_set], where archaic_set specifies 
which archaic individual(s) has/have the allele for it to be considered as "archaic". Options are 
"d_only" (found only in Denisovan), "n_only" (found in any of the Neanderthals but not the Denisovan),  
and "arch_all" (any archaic allele present in any individual).

This script generates a csv file, currently named "[archaic_set]_png_1k_venn.csv" (line 41). Each row
represents a pattern of archaic allele sharing, represented by a four digit binary value. The numbers 
represent presence or absence of the archaic allele in East Asians, Europeans, South Asians, and Papuans
(in order). A 0 means the allele is not found at least 1% frequency in a population from that superpopulation
(according to the 1000 Genomes nomenclature), while a 1 means it is present with at least 1% frequency. For 
example, "1111" means that the allele is found in all regions, while "0110" means that the allele is present
in Europeans and South Asians, but not East Asians or Papuans. The second column of the file is the count
of SNPs that show a given pattern. These numbers can be put into a program like venn_euler in R to make
a Venn Diagram.
"""
import gzip
import csv
import sys
import random
import re

archaic_set = sys.argv[1] #possible options: arch_all, n_only, d_only

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
vennOutfile = archaic_set + "_png_1k_venn.csv"

archLineList = ["aLine", "cLine", "vLine", "dLine"]
archLines = {"aLine":"", "cLine":"", "vLine":"", "dLine":""}
archPos = [0,0,0,0]

archLims = {"1": 249239825, "2": 243177959, "3":197874834, "4":191043233, "5":180715571, "6": 170919480, "7": 159128533, "8": 146300828, "9": 141102551, \
"10": 135512040, "11": 134946370, "12": 133841510, "13": 115108863, "14": 107289205, "15": 102506248, "16": 90174509, "17": 81160199,  \
"18": 78017036, "19": 59118836, "20": 62917942, "21": 48101478, "22": 51234514, "X": 155259633}
paupLims = {"1": 249240605, "2": 243188943, "3": 197961600, "4": 191044274, "5": 180904762, "6": 171053334, "7": 159128654, "8": 146304015, \
"9": 141149247, "10": 135524743, "11": 134946507, "12": 133841592, "13": 115108993, "14": 107289436, "15": 102521365, "16": 90294746, \
"17": 81195163, "18": 78017073, "19": 59118906, "20": 62965513, "21": 48119886, "22": 51244550, "X": 155260481}

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU", "PNG"]
afrPops = ["YRI", "LWK", "GWD", "MSL", "ESN"]

Pop_Tracker = {}
Archaic_Pattern_Venn = {}

for pop in populations:
    Pop_Tracker[pop] = []

for i in range(9,24):
    Pop_Tracker["PNG"].append(i)

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

min_sample_size = 15

for pop in Pop_Tracker.keys():
    if pop != "AFR" and pop != "PNG": #Papuans already only have 15 individuals
        ind_downsample = random.sample(Pop_Tracker[pop], min_sample_size)
        ind_downsample.sort()
        Pop_Tracker[pop] = ind_downsample

def decode_split_line(line):
    #line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
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
    non_afr_allele = "null"
    afr_alt_af = calc_afr_allele_freq(line)
    if afr_alt_af < 0.01 or afr_alt_af > 0.99:
        if afr_alt_af < 0.01: #allele of interest is alternate
            non_afr_allele = "1"
        else: #allele of interest is reference
            non_afr_allele = "0"
    return non_afr_allele

def arch_passes_quality_check(arch_info):
    arch_data_parse = arch_info.split(":")
    genotype_quality = arch_data_parse[7]
    if float(genotype_quality) >= 40.0:
        return True
    else:
        return False

def skip_archaic_header(aFile,line):
    while '#' in line:
        line = aFile.readline()
    return line

def scan_archaic_lines(aFile,archaicPos,afrPos):
    if chromosome in archLims:
        while int(archaicPos)<int(afrPos) and int(archaicPos)<archLims[chromosome]:
            line = aFile.readline()
            archspline = line.split()
            archaicPos = archspline[1]
    else:
        while int(archaicPos)<int(afrPos):
            line = aFile.readline()
            archspline = line.split()
            archaicPos = archspline[1]
    return line

def calc_pop_freq(population, allele):
    pop_allele_count = 0
    allele_total = len(Pop_Tracker[population]) * 2
    if chr == "X":
        allele_total = 0
    for ind in Pop_Tracker[population]:
        if population == "PNG":
            pop_allele_count += pSpline[ind].count(allele)
        else:
            pop_allele_count += mspline[ind].count(allele)
        if chr == "X": #count number of X alleles as we go to account for recombinant and non-recombinant
            if population == "PNG":
                allele_total = calc_total_X_allele(allele_total, pSpline[ind])
            else:
                allele_total = calc_total_X_allele(allele_total, mspline[ind])
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return rounded_freq

def correctArchSet(allele, genotypes): 
    fitsCriteria = False
    anea, cnea, vnea, den = genotypes[0:4]
    if archaic_set == "arch_all":
        fitsCriteria = any(allele in geno for geno in genotypes)
    elif archaic_set in ["d_only", "n_only"]:
        if all(geno != "-" for geno in genotypes):
            if archaic_set == "d_only":
                if allele in den and all(allele not in geno for geno in [anea, cnea, vnea]):
                    fitsCriteria = True
            elif archaic_set == "n_only":
                if allele not in den and any(allele in geno for geno in [anea, cnea, vnea]):
                    fitsCriteria = True        
    else:
        print("You have selected a non-standard archaic set")
    return fitsCriteria
      
for chromosome in chr_list:
    print(str(chromosome))
    kgFile = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    if chromosome == "X":
        kgFile = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    aNeaInfile = "/users/kwittdil/data/data/archaic_genomes/altai_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    cNeaInfile = "/users/kwittdil/data/data/archaic_genomes/chagyrskaya_VCF/ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/chr" + chromosome + ".noRB.vcf.gz"
    vNeaInfile = "/users/kwittdil/data/data/archaic_genomes/vindija_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    denInfile = "/users/kwittdil/data/data/archaic_genomes/denisova_VCF/ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr" + chromosome + "_mq25_mapab100.vcf.gz"
    papuanInfile = "./SGDP_Papuans_merged_chr" + str(chromosome) + ".vcf.gz"
    af = gzip.open(aNeaInfile, 'rt')
    cf = gzip.open(cNeaInfile, 'rt')
    vf = gzip.open(vNeaInfile, 'rt')
    df = gzip.open(denInfile, 'rt')
    pf = gzip.open(papuanInfile, 'rt')
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
    paupLine = pf.readline()
    paupLine = skip_archaic_header(pf,paupLine)
    pSpline = paupLine.split(sep="\t")
    paupPosition = pSpline[1]
    with gzip.open(kgFile, 'rt') as f:
        for line in f:
            mspline = decode_split_line(line)
            if '#' not in line:
                mod_position = mspline[1]
                while int(paupPosition) < int (mod_position) and int(paupPosition) != paupLims[chromosome]: #If there's a Papuan SNP not found in 1KG:
                    paupAlt = pSpline[4]
                    alternateAlleles = {paupAlt}
                    inArchaic = False
                    if int(paupPosition) <= archLims[chromosome]:
                        arch_genotypes = ["-","-","-","-"]
                        for j in range(0, len(archFiles)):
                            archaicPosition = archPos[j]
                            archaicLine = archLineList[j]
                            currArchLine = archLines[archaicLine]
                            afile = archFiles[j]
                            if int(archaicPosition) < int(paupPosition):
                                archLines[archaicLine] = currArchLine = scan_archaic_lines(afile,archaicPosition,paupPosition)
                            aSpline = currArchLine.split()
                            archaicPosition,_,arch_ref,arch_alt,genotype_qual = aSpline[1:6]
                            if arch_alt != ".":
                                alternateAlleles.add(arch_alt)
                            archPos[j] = archaicPosition
                            if int(archaicPosition) == int(paupPosition):
                                inArchaic = True
                                if float(genotype_qual) >= 40:
                                    arch_genotypes[j] = aSpline[9][0:3]
                        if len(alternateAlleles) == 1:
                            derivedFreq = calc_pop_freq("PNG", "1")
                            if correctArchSet("1", arch_genotypes) and inArchaic: # If archaic humans also have the allele
                                if derivedFreq > 0:
                                    vennPattern = "0001"
                                    if vennPattern not in Archaic_Pattern_Venn:
                                        Archaic_Pattern_Venn[vennPattern] = 1
                                    else:
                                        Archaic_Pattern_Venn[vennPattern] += 1
                    paupLine = pf.readline()
                    pSpline = paupLine.split()
                    paupPosition = pSpline[1]
                if snp_is_biallelic(mspline):
                    mod_alt_allele = mspline[4]
                    alternateAlleles = {mod_alt_allele}
                    if int(paupPosition) == int(mod_position):
                        paupAlt = pSpline[4]
                        if paupAlt != ".":
                            alternateAlleles.add(paupAlt)
                    non_afr_allele = identify_non_afr_alleles(mspline)  
                    if non_afr_allele in ["0", "1"]:  
                        freq_by_pop = []
                        arch_genotypes = ["-","-","-","-"]
                        inArchaic = False
                        if int(mod_position) <= archLims[chromosome]:
                            for j in range(0, len(archFiles)):
                                archaicPosition = archPos[j]
                                archaicLine = archLineList[j]
                                currArchLine = archLines[archaicLine]
                                afile = archFiles[j]
                                if int(archaicPosition) < int(mod_position):
                                    archLines[archaicLine] = currArchLine = scan_archaic_lines(afile,archaicPosition,mod_position)
                                aSpline = currArchLine.split()
                                archaicPosition,_,arch_ref,arch_alt,genotype_qual = aSpline[1:6]
                                if arch_alt != ".":
                                    alternateAlleles.add(arch_alt)
                                archPos[j] = archaicPosition
                                if int(archaicPosition) == int(mod_position):
                                    if float(genotype_qual) >= 40:
                                        arch_genotypes[j] = aSpline[9][0:3]
                                        inArchaic = True
                            if len(alternateAlleles) == 1:
                                vennPattern = ""
                                for pop in populations[:-1]:
                                    nonafr_freq = calc_pop_freq(pop, non_afr_allele)
                                    freq_by_pop.append(nonafr_freq)
                                if correctArchSet(non_afr_allele,arch_genotypes) and inArchaic: 
                                    vennPattern = vennPattern + ("1" if any(x>0.01 for x in freq_by_pop[1:6]) else "0")  #EAS
                                    vennPattern = vennPattern + ("1" if any(x>0.01 for x in freq_by_pop[6:11]) else "0") #EUR
                                    vennPattern = vennPattern + ("1" if any(x>0.01 for x in freq_by_pop[15:20]) else "0") #SAS
                                    if int(paupPosition) == int(mod_position):
                                        derivedFreq = calc_pop_freq("PNG", "1")
                                        if non_afr_allele == "1":
                                            nonafr_freq = derivedFreq
                                        else:
                                            nonafr_freq = 1-derivedFreq
                                        vennPattern = vennPattern + "1" if nonafr_freq>0.01 else vennPattern + "0" #PNG
                                    else:
                                        vennPattern += "0"
                                    if vennPattern not in Archaic_Pattern_Venn:
                                        Archaic_Pattern_Venn[vennPattern] = 1
                                    else:
                                        Archaic_Pattern_Venn[vennPattern] += 1
                if int(paupPosition) == int(mod_position) and int(paupPosition) < paupLims[chromosome]:
                    paupLine = pf.readline()
                    pSpline = paupLine.split()
                    paupPosition = pSpline[1]
    while int(paupPosition) != paupLims[chromosome]: #If there's a Papuan SNP not found in 1KG:
        paupAlt = pSpline[4]
        alternateAlleles = {paupAlt}
        if int(paupPosition) <= archLims[chromosome]:
            arch_genotypes = ["-","-","-","-"]
            inArchaic = False
            for j in range(0, len(archFiles)):
                archaicPosition = archPos[j]
                archaicLine = archLineList[j]
                currArchLine = archLines[archaicLine]
                afile = archFiles[j]
                if int(archaicPosition) < int(paupPosition):
                    archLines[archaicLine] = currArchLine = scan_archaic_lines(afile,archaicPosition,paupPosition)
                aSpline = currArchLine.split()
                archaicPosition,_,arch_ref,arch_alt,genotype_qual = aSpline[1:6]
                alternateAlleles.add(arch_alt)
                archPos[j] = archaicPosition
                if int(archaicPosition) == int(paupPosition):
                    if float(genotype_qual) >= 40:
                        inArchaic = True
                        arch_genotypes[j] = aSpline[9][0:3]
            if len(alternateAlleles) == 1:
                derivedFreq = calc_pop_freq("PNG", "1")
                if correctArchSet("1", arch_genotypes) and inArchaic: # If archaic humans also have the allele
                    if derivedFreq > 0:
                        vennPattern = "0001"
                        if vennPattern not in Archaic_Pattern_Venn:
                            Archaic_Pattern_Venn[vennPattern] = 1
                        else:
                            Archaic_Pattern_Venn[vennPattern] += 1
        paupLine = pf.readline()
        pSpline = paupLine.split()
        paupPosition = pSpline[1]    

with open(vennOutfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["pattern", "count"]
    w.writerow(header)
    for pattern in Archaic_Pattern_Venn:
        allele_line = [pattern,str(Archaic_Pattern_Venn[pattern])]
        w.writerow(allele_line)