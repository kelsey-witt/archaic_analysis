"""
This script was written by Kelsey Witt Dillon in April 2021. The goal of this script is to 
make a file containing all biallelic sites from the 1000 Genomes dataset where one allele is
at less than 1% frequency in Africans, at least 1% frequency in at least one population outside
of Africa, and shared with sequenced archaic genomes (Altai Chagyrskaya and Vindija Neanderthal
and Denisovan). It records whether or not each non-African individual has the putatively archaic
allele in that position.

The input files are the csv files for each chromosome containing alleles that are rare in Africa 
(which come from ID_rare_AFR_snps.py, "afrInfile", line 62), and the four archaic human VCF
files ("aNeaInfile", "vNeaInfile","cNeaInfile", "denInfile", lines 63-66.)

Usage: python3 ID_rare_AFR_archaic_snps.py

This script generates a csv file, currently named "archaic_rareAFR_snps.csv" (line 27). Each row represents a 
site in the genome, and the first nine columns are the chromosome number, position, the reference allele, 
the alternate allele, and the allele that is rare in Africa (0 if reference, 1 if alternate), and the
genotypes of the Altai Neanderthal, Chagyrskaya Neanderthal, Vindija Neanderthal, and Denisovan. The remaining 
columns are the presence or absence of the possibly archaic allele in that individual's genotype, organized 
by population. A "0" means that the individual lacks the archaic allele, while a "1" indicates that the 
archaic allele is present in that individual. The header lists the individual's ID as their population 
and the 0-indexed column number they appear in the 1000 Genomes vcf, separated by an underscore (ie "GBR_9").
"""
import gzip
import csv

outfile = "archaic_rareAFR_snps.csv"

archLineList = ["aLine", "cLine", "vLine", "dLine"]
archLines = {"aLine":"", "cLine":"", "vLine":"", "dLine":""}
archPos = [0,0,0,0]
archLims = {"2": 243177959, "3":197874834, "4":191043233, "5":180715571, "6": 170919480, "7": 159128533, "8": 146300828, "9": 141102551, \
    "10": 135512040, "11": 134946370, "12": 133841510, "13": 115108863, "14": 107289205, "15": 102506248, "16": 90174509, "17": 81160199,  \
        "18": 78017036, "20": 62917942, "21": 48101478, "22": 51234514}

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

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))
chr_list.append("X")

with open(outfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")         
    for chromosome in chr_list:
        afrInfile = "rare_AFR_snps_chr" + chromosome + ".csv.gz"
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
                if "Chromosome" in line:
                    header = ["Chromosome","Position","Reference","Alternate","Rare Allele","Altai","Chagyrskaya","Vindija","Denisovan"]
                    spline = line[:-1].split(sep=",")
                    for colnum in range(5,len(spline)):
                        header.append(spline[colnum])
                    w.writerow(header)
                else:
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
                        rareInArch = False
                        for genotype in arch_genotypes:
                            if rare_allele in genotype:
                                rareInArch = True
                        if rareInArch:
                            outLine = spline[0:5] + arch_genotypes + spline[5:]
                            w.writerow(outLine)