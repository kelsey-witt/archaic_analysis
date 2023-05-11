"""
This script was written by Kelsey Witt Dillon in December 2021. The goal of this script is to 
count the number of archaic alleles found in each non-African individual in the 1000 Genomes
populations and the Papuans in the Simons Genome Diversity Project, and to also calculate 
allele frequencies for both archaic and modern (not shared with archaic genomes) alleles with 
low frequency (<1%) in Africa, and classify them as rare (above 1% frequency but below 20% 
frequency) or common (20% frequency or greater). For an allele to be classified as archaic, 
it had to be found in Africans at a frequency less than 1% but at a frequency of 1% or more 
in at least one non-African population and shared with the archaic humans of interest (as 
established by archaic_set).

The input files are the .panel file and VCF files from the 1000 Genomes Phase 3 dataset (pop_file 
in line 47 and kgFile in line 243), a merged vcf containing the genotype calls for the SGDP Papuans
separated by chromosome (papuanInfile, line 248) and the four archaic human VCF files ("aNeaInfile", 
"vNeaInfile","cNeaInfile", "denInfile", lines 244-247.)

Usage: python3 archaic_allele_counter_mp_papuans_n15.py [archaic_set], where archaic_set specifies 
which archaic individual(s) has/have the allele for it to be considered as "archaic". Options are 
"d_only" (found only in Denisovan), "n_only" (found in any of the Neanderthals but not the Denisovan),  
and "arch_all" (any archaic allele present in any individual).

This script generates two output files. The first, "[archaic_set]_png_1k_ctable.csv" (popOutfile, 
line 45) is a csv file that summarizes population-level information. Each row is a population, 
and the columns are, in order, the population name, the number of non-archaic alleles that are rare in 
Africa but 20%+ frequency in that population ("modern common"), the number of non-archaic alleles that 
are rare in Africa but between 1 and 20% frequency in that population ("modern rare"), the number of 
archaic alleles that are rare in Africa but 20%+ frequency in that population ("archaic common"), and 
the number of archaic alleles that are rare in Africa but between 1 and 20% frequency in that population 
("archaic rare"). 

The second output file, "[archaic_set]_png_1k_archaic_score.csv" (indOutfile, line 46) is a 
csv file that summarizes individual-level information. Each row represents an individual, and the 
columns are, in order, the individual ID, the population the individual is from, and the number of 
archaic alleles that individual has. "Archaic" alleles are determined by the criteria above, and are 
counted across all sites in the genome (ie a site that is homozygous for an archaic allele would 
increase the total by 2).
"""
import gzip
import csv
import sys
import random
import re

archaic_set = sys.argv[1] #possible options: arch_all, n_only, d_only

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
popOutfile = archaic_set + "_png_1k_ctable.csv"
indOutfile = archaic_set + "_png_1k_archaic_score.csv"

#pop_file="./toy_popfile.txt"
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
#chr_list.append("X")

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU", "PNG"]
afrPops = ["YRI", "LWK", "GWD", "MSL", "ESN"]

Pop_Tracker = {}
Male_Pop_Tracker = {}
NonAfr_Archaic_Allele_Counter = {}
NonAfr_Modern_Allele_Counter = {}
Ind_Tracker = {}
Ind_Allele_Counter = {}
Seg_Site_counter = {}

for pop in populations:
    Pop_Tracker[pop] = []
    if pop != "AFR":
        NonAfr_Modern_Allele_Counter[pop] = {}
        NonAfr_Modern_Allele_Counter[pop]["common"] = 0
        NonAfr_Modern_Allele_Counter[pop]["rare"] = 0
        NonAfr_Archaic_Allele_Counter[pop] = {}
        NonAfr_Archaic_Allele_Counter[pop]["common"] = 0
        NonAfr_Archaic_Allele_Counter[pop]["rare"] = 0
        Seg_Site_counter[pop] = 0

for i in range(9,24):
    Pop_Tracker["PNG"].append(i)
    Ind_Allele_Counter[i+5000] = 0

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
            Ind_Tracker[col_counter] = line_col[0]
            Ind_Allele_Counter[col_counter] = 0
        col_counter += 1

Ind_Tracker.update({5009:"LP6005441-DNA_A10", 5010:"LP6005441-DNA_B10", 5011:"LP6005443-DNA_A08", 5012:"LP6005443-DNA_B08", 5013:"LP6005443-DNA_C07", \
5014:"LP6005443-DNA_C08", 5015:"LP6005443-DNA_D07", 5016:"LP6005443-DNA_D08", 5017:"LP6005443-DNA_E07", 5018:"LP6005443-DNA_E08", 5019:"LP6005443-DNA_F07",\
5020:"LP6005443-DNA_F08", 5021:"LP6005443-DNA_G07", 5022:"LP6005443-DNA_H07", 5023:"SS6004472"})

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

def allele_counter(freq_list, counter_dict):
    if freq_list[0] < 0.01:
        for x in range(1,len(freq_list)):
            if float(freq_list[x]) > 0.20:
                counter_dict[populations[x]]["common"] += 1
            else:
            #elif float(freq_list[x]) > 0.01:
                counter_dict[populations[x]]["rare"] += 1

def allele_counter_papuan(freq, counter_dict):
    if float(freq) > 0.20:
        counter_dict["PNG"]["common"] += 1
    elif float(freq) > 0.04:
        counter_dict["PNG"]["rare"] += 1

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
                    alleleCounter = 0
                    allele_total = 32
                    if chromosome == "X":
                        allele_total = 0
                        for ind in Pop_Tracker["PNG"]:
                            allele_total = calc_total_X_allele(allele_total, pSpline[ind])
                    for ind in Pop_Tracker["PNG"]:
                        alleleCounter += pSpline[ind][0:3].count("1")
                    if alleleCounter > 0 and alleleCounter < allele_total:
                        Seg_Site_counter["PNG"] += 1
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
                                allele_counter_papuan(derivedFreq, NonAfr_Archaic_Allele_Counter)
                                for ind in Pop_Tracker["PNG"]:
                                    Ind_Allele_Counter[ind+5000] += pSpline[ind][0:3].count("1")
                            else:
                                allele_counter_papuan(derivedFreq, NonAfr_Modern_Allele_Counter)
                    paupLine = pf.readline()
                    pSpline = paupLine.split()
                    paupPosition = pSpline[1]
                if snp_is_biallelic(mspline):
                    mod_alt_allele = mspline[4]
                    for pop in populations[1:-1]:
                        alleleCounter = 0
                        allele_total = len(Pop_Tracker[pop]) * 2
                        if chromosome == "X":
                            allele_total=0
                            for ind in Pop_Tracker[pop]:
                                allele_total = calc_total_X_allele(allele_total, mspline[ind])
                        for ind in Pop_Tracker[pop]:
                            alleleCounter += mspline[ind].count("1")
                        if alleleCounter > 0 and alleleCounter < allele_total: # as long as they aren't all 0s or 1s it's a segregating site
                            Seg_Site_counter[pop] += 1
                    alternateAlleles = {mod_alt_allele}
                    if int(paupPosition) == int(mod_position):
                        paupAlt = pSpline[4]
                        if paupAlt != ".":
                            alternateAlleles.add(paupAlt)
                        alleleCounter = 0
                        allele_total = 32
                        if chromosome == "X":
                            allele_total=0
                            for ind in Pop_Tracker["PNG"]:
                                allele_total = calc_total_X_allele(allele_total, pSpline[ind])
                        for ind in Pop_Tracker["PNG"]:
                            alleleCounter += pSpline[ind][0:3].count("1")
                        if alleleCounter > 0 and alleleCounter < allele_total: # as long as they aren't all 0s or 1s it's a segregating site
                            Seg_Site_counter["PNG"] += 1
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
                                for pop in populations[:-1]:
                                    nonafr_freq = calc_pop_freq(pop, non_afr_allele)
                                    freq_by_pop.append(nonafr_freq)
                                if correctArchSet(non_afr_allele,arch_genotypes) and inArchaic: 
                                    allele_counter(freq_by_pop, NonAfr_Archaic_Allele_Counter)
                                    for ind in Pop_Tracker[pop]:
                                            Ind_Allele_Counter[ind] += mspline[ind].count(non_afr_allele)
                                else:
                                    allele_counter(freq_by_pop, NonAfr_Modern_Allele_Counter)
                                if int(paupPosition) == int(mod_position):
                                    derivedFreq = calc_pop_freq("PNG", "1")
                                    if non_afr_allele == "1":
                                        nonafr_freq = derivedFreq
                                    else:
                                        nonafr_freq = 1-derivedFreq
                                    if correctArchSet(non_afr_allele, arch_genotypes) and inArchaic: # If archaic humans also have the allele
                                        allele_counter_papuan(nonafr_freq, NonAfr_Archaic_Allele_Counter)
                                        for ind in Pop_Tracker["PNG"]:
                                            Ind_Allele_Counter[ind+5000] += pSpline[ind][0:3].count("1")
                                    else:
                                        allele_counter_papuan(nonafr_freq, NonAfr_Modern_Allele_Counter)
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
                    allele_counter_papuan(derivedFreq, NonAfr_Archaic_Allele_Counter)
                    for ind in Pop_Tracker["PNG"]:
                        Ind_Allele_Counter[ind+5000] += pSpline[ind][0:3].count("1")
                else:
                    allele_counter_papuan(derivedFreq, NonAfr_Modern_Allele_Counter)
        paupLine = pf.readline()
        pSpline = paupLine.split()
        paupPosition = pSpline[1]    
    
with open(indOutfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["ind", "pop", "positions"]
    w.writerow(header)
    for pop in populations[1:-1]:
        for ind in Pop_Tracker[pop]:
            indID = Ind_Tracker[ind]
            indScore = str(Ind_Allele_Counter[ind])
            score_line = [indID,pop,indScore]
            w.writerow(score_line)
    for ind in Pop_Tracker["PNG"]:
        higherVal = 5000 + ind
        indID = Ind_Tracker[higherVal]
        indScore = str(Ind_Allele_Counter[higherVal])
        score_line = [indID,"PNG",indScore]
        w.writerow(score_line)

with open(popOutfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["pop", "mod_common", "mod_rare", "archaic_common", "archaic_rare", "seg_sites"]
    w.writerow(header)
    for i in range(1,len(populations)):
        allele_line = [populations[i], str(NonAfr_Modern_Allele_Counter[populations[i]]["common"]), str(NonAfr_Modern_Allele_Counter[populations[i]]["rare"]),  str(NonAfr_Archaic_Allele_Counter[populations[i]]["common"]), str(NonAfr_Archaic_Allele_Counter[populations[i]]["rare"]), str(Seg_Site_counter[populations[i]])]
        w.writerow(allele_line)