"""
This script was written by Kelsey Witt Dillon in October 2021. The goal of this script is to 
take the archaic ancestry tracts published by Sankararaman et al. 2016 and assess how archaic
ancestry coverage in East Asians vs. Europeans changes with increasing sample size. The script 
first combines the haploid ancestry tracts from each individual into diploid ancestry tracts
(where archaic ancestry on either chromosome counts as diploid archaic ancestry) and then
sums the total length of tracts with archaic ancestry for a sample of individuals. Sample sizes 
range from 1-20 (and any regions with more than 20 individuals are downsized to that number), 
and the script samples 100 replicates for each sample size.

The input files are the .ind and .haplotypes files from the Sankararaman et al. 2016 ancestry
tracts (indFile, line 48 and hapFile, line 68)

Usage: python3 sank_interval_summing.py

This script generates a csv file labeled "sank_map_eur_eas_cvtable.csv", which is a table of
the amount of Neanderthal and Denisovan coverage in 100 replicates of each sample size for each 
geographic region. Each row represents a replicate of one region/archaic source/sample size 
combination, and the columns are, in order, the replicate number for that combination, the 
sample size, the geographic region (eastasia or westeurasia as defined by SGDP), the archaic 
source(neanderthal or denisovan), and the total length of genomic regions for that sample 
of individuals that have archaic ancestry.
"""

import random

outfile = "./sank_map_eur_eas_cvtable.csv"

h=open(outfile,'w')

sample_sizes = [1,5,10,20]
archaic_pops = ["denisova","neandertal"]
regions = ["eastasia","westeurasia"] #can be any combination of regions from the Simons Genome Diversity Project data
chromosomes = ["X"]
for i in range(1,23):
    chromosomes.append(str(i))

outHeader = ["rep","sample_size","region","archaic","tract_len"]
outLine = ",".join(outHeader) + "\n"
h.write(outLine)

for region in regions:
    indFile = "./sank_introgression_maps/sgdp/ids/" + region + ".ind"
    with open(indFile) as f:
        lineCount = 0
        indList = []
        for line in indFile:
            if lineCount%2 ==0:
                indList.append(lineCount//2)
            lineCount += 1
    for ss in sample_sizes:
        for rep in range(0,100):
            indSubsample = random.sample(indList, ss)
            lineList = []
            for ind in indSubsample:
                lineList.append(ind*2)
                lineList.append(ind*2+1)
            for archaic in archaic_pops:
                tractLen = 0
                for chromosome in chromosomes:
                    popTracts = []
                    mergedTracts = []
                    hapFile = "./sank_introgression_maps/sgdp/2/" + archaic + "/" + region + "/summaries/haplotypes/" + chromosome + ".thresh-50.length-0.00.haplotypes"
                    with open(hapFile) as g:
                        for line in g:
                            if '#' not in line:
                                spline = line.split()
                                ind,start,end = spline[1:4]
                                if int(ind)//2 in indSubsample:
                                    popTracts.append([int(start),int(end)])
                    if popTracts:
                        tractSort = sorted(popTracts)
                        for i in range(0,len(tractSort)):
                            if i == 0:
                                tractSt = int(tractSort[i][0])
                                tractEnd = int(tractSort[i][1])
                            else:
                                currTractSt = int(tractSort[i][0])
                                currTractEnd = int(tractSort[i][1])
                                if tractEnd < currTractSt:
                                    mergedTracts.append([tractSt,tractEnd])
                                    tractSt = currTractSt
                                    tractEnd = currTractEnd
                                elif tractEnd < currTractEnd:
                                    tractEnd = currTractEnd
                        mergedTracts.append([tractSt,tractEnd])
                        for tract in mergedTracts:
                            trLen = tract[1] - tract[0]
                            tractLen += trLen
                outCols=[str(rep),str(ss),region,archaic,str(tractLen)]
                outLine = ",".join(outCols)+"\n"
                h.write(outLine)
h.close()
