import gzip
import re
import sys
import csv

archaic_set = sys.argv[1] #possible options: d_only, n_only, d_all, n_all, nd_both, nd_either

#pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
outfile = archaic_set + "_sprime_archaic_cvg.csv"

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Pop_Tracker = {}

for pop in populations:
    Pop_Tracker[pop] = []

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        pop = line_col[1]
        if pop in populations:
            pop = line_col[1]
            Pop_Tracker[pop].append(col_counter)
        col_counter += 1

def decode_split_line(line):
    line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

with open(outfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["Chromosome","Position","Reference","Alternate","Archaic Allele"]
    for pop in populations:
        for ind in Pop_Tracker[pop]:
            indTitle = pop + "_" + str(ind)
            header.append(indTitle)
    w.writerow(header)
    for chromosome in chr_list:
        #modern_file = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        modern_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        alleleMatch = {}
        allArchaicSNPs = set()
        for pop in populations:
            sprime_file = "./sprime_data/" + pop + ".chr" + chromosome + ".ND_match"
            #sprime_file = "./Sprime/perchrom_files/" + pop + ".chr" + chromosome + ".ND_match"
            alleleMatch[pop] = {}
            with open(sprime_file, 'r') as sf:
                next(sf)
                for line in sf:
                    line_split = line.split()
                    position = line_split[1]
                    arch_allele = line_split[6] 
                    Nstatus, Dstatus = line_split[8:10]
                    isArchaic = False
                    if Nstatus == "match" or Dstatus == "match":
                        if archaic_set == "arch_all":
                            isArchaic = True
                        elif archaic_set == "n_only":
                            if Nstatus == "match" and Dstatus == "mismatch":
                                isArchaic = True
                        elif archaic_set == "d_only":
                            if Dstatus == "match" and Nstatus == "mismatch":
                                isArchaic = True
                    if isArchaic:
                        allArchaicSNPs.add(position)
                        alleleMatch[pop][position] = arch_allele
        with gzip.open(modern_file, 'r') as mf:
            for mod_line in mf:
                mod_line = decode_split_line(mod_line)
                mod_position = mod_line[1]
                if mod_position in allArchaicSNPs:
                    ref_allele,alt_allele = mod_line[3:5]
                    outLine = [str(chromosome),str(position),ref_allele,alt_allele]
                    for pop in populations:
                        if mod_position in alleleMatch[pop].keys():
                            arch_allele = alleleMatch[pop][mod_position]
                            if len(outLine)==4:
                                outLine.append(str(arch_allele))
                    for pop in populations:
                        if mod_position in alleleMatch[pop].keys():
                            for ind in Pop_Tracker[pop]:
                                if arch_allele in mod_line[ind]:
                                    outLine.append("1")
                                else:
                                    outLine.append("0")
                        else:
                            for ind in Pop_Tracker[pop]:
                                outLine.append("0")
                    w.writerow(outLine)

