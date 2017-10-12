import gzip
import re 

#pop_file = "/home/1000genome_data/integrated_call_samples_v3.20130502.ALL.panel"
pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
outfile = "Nea_derived_ctable_pop.txt"

allele_present = ["0/1", "1/1"]
populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Pop_Tracker = {}
Male_Pop_Tracker = {}

for pop in populations:
    Pop_Tracker[pop] = []
    Male_Pop_Tracker[pop] = []

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))
chr_list.append("X")

allele_present = ["0/1", "1/1"]

col_counter = 9
male_col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] == "AFR":
            Pop_Tracker["AFR"].append(col_counter)
        else:
            pop = line_col[1]
            Pop_Tracker[pop].append(col_counter)
        col_counter += 1
        if line_col[3] == "male":
            if line_col[2] == "AFR":
                Male_Pop_Tracker["AFR"].append(male_col_counter)
            else:
                pop = line_col[1]
                Male_Pop_Tracker[pop].append(male_col_counter)
            male_col_counter += 1
print(Male_Pop_Tracker)
print(male_col_counter)

Archaic_Allele_Counter = {}
Modern_Allele_Counter = {}
for pop in populations: 
    Modern_Allele_Counter[pop] = {}
    Modern_Allele_Counter[pop]["common"] = 0
    Modern_Allele_Counter[pop]["rare"] = 0
    Archaic_Allele_Counter[pop] = {}
    Archaic_Allele_Counter[pop]["common"] = 0
    Archaic_Allele_Counter[pop]["rare"] = 0

"""
class Allele:
    def __init__(self, reference, alternate, ancestral):
        self.reference = reference
        self.alternate = alternate
        self.ancestral = ancestral
    
    def 

x = Allele('A', 'C', 'G')
# class AlleleCounters:
#     def __init__(self):
#         self.common = 0
#         self.rare = 0
"""

def pop_freq(population, derived_allele):
    pop_allele_count = 0
    for ind in Pop_Tracker[population]:
        num_alt = 0
        if len(mod_line_col[ind]) == 1:
            if mod_line_col[ind] == derived_allele:
                num_alt = 1
        else:
            if mod_line_col[ind] == "1|0" or mod_line_col[ind] == "0|1":
                num_alt = 1
            elif mod_line_col[ind] == "1|1" and derived_allele == 1:
                num_alt = 2
            elif mod_line_col[ind] =="0|0" and derived_allele == 0:
                num_alt = 2
        pop_allele_count += num_alt
    allele_total = len(Pop_Tracker[population]) * 2
    derived_frequency = float(pop_allele_count/allele_total)
    return derived_frequency

def pop_freq_y(population, derived_allele):
    pop_allele_count = 0
    for ind in Male_Pop_Tracker[population]:
            num_alt = 0
            if mod_line_col[ind] == derived_allele:
                num_alt = 1
            pop_allele_count += num_alt
    allele_total = len(Male_Pop_Tracker[population]) * 2
    derived_frequency = float(pop_allele_count/allele_total)
    return derived_frequency


def allele_counter(freq_list, counter_dict):
    if freq_list[0] < 0.01: #if it's not found in Africans
        for x in range(1,len(freq_list)):
            if freq_list[x] > 0.20:
                counter_dict[populations[x]]["common"] += 1
            elif freq_list[x] > 0.0:
                counter_dict[populations[x]]["rare"] += 1

for d in range(0, len(chr_list)):
    modern_infile = "/home/1000genome_data/ALL.chr" + chr_list[d] + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    archaic_infile = "./Nea_Den_" + chr_list[d] + "_derived.txt.gz"
    Neanderthal_Variation = {}
    with gzip.open(archaic_infile) as a:
        for line in a:
            line = str(line.decode('utf-8'))
            arch_line_col = line.split()
            if arch_line_col[4] in allele_present:
                position, _, arch_allele, chimp_allele = arch_line_col[0:4]
                if arch_allele != chimp_allele and arch_allele != ".":
                    Neanderthal_Variation[position] = chimp_allele
    with gzip.open(modern_infile) as b:
        for line in b:
            line = str(line.decode('utf-8'))
            derived_freq_by_pop = []
            mod_line_col = line.split()
            # use itertools.islice to skip header lines
            if "#" not in line and mod_line_col[1] in Neanderthal_Variation.keys():
                if Neanderthal_Variation[mod_line_col[1]] == mod_line_col[3]:
                    for pop in populations:
                        der_freq = pop_freq(pop,1)
                        derived_freq_by_pop.append(der_freq)
                    allele_counter(derived_freq_by_pop, Archaic_Allele_Counter)
                else:
                    for pop in populations:
                        der_freq = pop_freq(pop,0)
                        derived_freq_by_pop.append(der_freq)
                    allele_counter(derived_freq_by_pop, Archaic_Allele_Counter)
            elif "#" not in line and len(mod_line_col[3]) == 1 and len(mod_line_col[4]) == 1:
                if re.match(".+;AA=(\S)|", line) is not None:
                    ancestor = re.search(".+;AA=(\S)|", mod_line_col[7])
                    ancestral_allele = ancestor.group(1)
                    if ancestral_allele == ".":
                        for pop in populations:
                            der_freq = pop_freq(pop,1)
                            derived_freq_by_pop.append(der_freq)
                        allele_counter(derived_freq_by_pop, Modern_Allele_Counter)
                    elif ancestral_allele == mod_line_col[4]:
                        for pop in populations:
                            der_freq = pop_freq(pop,0)
                            derived_freq_by_pop.append(der_freq)
                        allele_counter(derived_freq_by_pop, Modern_Allele_Counter)

modern_infile = "/home/1000genome_data/ALL.chrY.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
archaic_infile = "./Nea_Den_Y_derived.txt.gz"
Neanderthal_Variation = {}
with gzip.open(archaic_infile) as f:
    for line in f:
        line = str(line.decode('utf-8'))
        arch_line_col = line.split()
        if arch_line_col[4] in allele_present:
            position, _, arch_allele, chimp_allele = arch_line_col[0:4]
            if arch_allele != chimp_allele and arch_allele != ".":
                Neanderthal_Variation[position] = chimp_allele
with gzip.open(modern_infile) as g:
    for line in g:
        line = str(line.decode('utf-8'))
        derived_freq_by_pop = []
        mod_line_col = line.split()
            # use itertools.islice to skip header lines
        if "#" not in line and mod_line_col[1] in Neanderthal_Variation.keys():
            if Neanderthal_Variation[mod_line_col[1]] == mod_line_col[3]:
                for pop in populations:
                    der_freq = pop_freq_y(pop, 1)
                    derived_freq_by_pop.append(der_freq)
                allele_counter(derived_freq_by_pop, Archaic_Allele_Counter)
            else:
                for pop in populations:
                    der_freq = pop_freq_y(pop, 1)
                    derived_freq_by_pop.append(der_freq)
                allele_counter(derived_freq_by_pop, Archaic_Allele_Counter)
        elif "#" not in line and len(mod_line_col[3]) == 1 and len(mod_line_col[4]) == 1:
            if re.match(".+;AA=(\S)|", line) is not None:
                ancestor = re.search(".+;AA=(\S)|", mod_line_col[7])
                ancestral_allele = ancestor.group(1)
                if ancestral_allele == ".":
                    for pop in populations:
                        der_freq = pop_freq_y(pop,1)
                        derived_freq_by_pop.append(der_freq)
                    allele_counter(derived_freq_by_pop, Modern_Allele_Counter)
                elif ancestral_allele == mod_line_col[4]:
                    for pop in populations:
                        der_freq = pop_freq_y(pop,0)
                        derived_freq_by_pop.append(der_freq)
                    allele_counter(derived_freq_by_pop, Modern_Allele_Counter)


with open(outfile, 'w') as text_file:
    population_headers = ["African Populations", "Chinese", "Japanese", "Southern Han Chinese", "Chinese Dai", "Vietnamese", "Western Europeans", "Italians", "Finnish", "British", "Spanish", "Mexicans", "Puerto Ricans", "Colombians", "Peruvians", "Gujarati", "Punjabi", "Bengali", "Sri Lankans", "Telugu"]
    for i in range(1,len(populations)):
        text_file.write(population_headers[i] + "\n")
        text_file.write("Common Modern Alleles: " + str(Modern_Allele_Counter[super_populations[i]]["common"]) + "\n")
        text_file.write("Rare Modern Alleles: " + str(Modern_Allele_Counter[super_populations[i]]["rare"]) + "\n")
        text_file.write("Common Archaic Alleles: " + str(Archaic_Allele_Counter[super_populations[i]]["common"]) + "\n")
        text_file.write("Rare Archaic Alleles: " + str(Archaic_Allele_Counter[super_populations[i]]["rare"]) + "\n")
        if Archaic_Allele_Counter[super_populations[i]]["common"] !=0:
            text_file.write("Archaic rare:common ratio: " + str(float(Archaic_Allele_Counter[super_populations[i]]["rare"]/Archaic_Allele_Counter[super_populations[i]]["common"])) + "\n")
        if Modern_Allele_Counter[super_populations[i]]["common"] !=0:
            text_file.write("Modern rare:common ratio: " + str(float(Modern_Allele_Counter[super_populations[i]]["rare"]/Modern_Allele_Counter[super_populations[i]]["common"])) + "\n")
        if Archaic_Allele_Counter[super_populations[i]]["common"] !=0 and Modern_Allele_Counter[super_populations[i]]["common"] != 0:
            archaic_freq = float(Archaic_Allele_Counter[super_populations[i]]["rare"]/Archaic_Allele_Counter[super_populations[i]]["common"])
            modern_freq = float(Modern_Allele_Counter[super_populations[i]]["rare"]/Modern_Allele_Counter[super_populations[i]]["common"])
            modern_archaic_ratio = archaic_freq/modern_freq
            text_file.write("Archaic and Modern Frequency Ratio: " + str(modern_archaic_ratio) + "\n \n")
        else:
            text_file.write("Could not calculate ratio - division by zero \n \n")