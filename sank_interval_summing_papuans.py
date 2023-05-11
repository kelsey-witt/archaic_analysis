"""
This script was written by Kelsey Witt Dillon in January 2021. The goal of this script is to 
identify the combined length of archaic introgression tracts in a subset of individuals from 
populations from the Simons Genome Diversity Project, organized by geographic region. Haploid 
archaic ancestry tracts are combined and summed over the sample of individuals for a total
length of tracts with archaic ancestry. Sample sizes range from 1-16 (and any regions with 
more than 16 individuals are downsized to that number), and the script samples 100 replicates 
for each sample size. 

The input files are the individual ID files and archaic haplotype files from Sankararaman et al.
2016 (indFile in line 45 and hapFile in line 69)

Usage: python3 sank_interval_summing_papuans.py

This script generates a csv file labeled "sank_map_png_cvtable_papuans.csv", which is a table of
the amount of Neanderthal and Denisovan coverage in 100 replicates of each sample size for each 
geographic region. Each row represents a replicate of one region/archaic source/sample size combination,
and the columns are, in order, the replicate number for that combination, the sample size,
the geographic region (eastasia, oceania, southasia, or westeurasia as defined by SGDP), the
archaic source(neanderthal or denisovan), and the total length of genomic regions for that sample of 
individuals that have archaic ancestry.
"""

import random

outfile = "./sank_map_png_cvtable_papuans.csv"

h=open(outfile,'w')

sample_sizes = [1,5,10,16]
archaic_pops = ["denisova","neandertal"]
regions = ["eastasia","oceania","southasia","westeurasia"]
pops = {"eastasia": ["Dai","Daur","Han","Hezhen","Lahu","Miao","Naxi","Oroqen","She","Tu","Tujia","Uygur","Xibo","Yi"], "oceania":["Papuan"], \
"southasia":["Brahmin","Irula","Bengali","Kapu","Madiga","Mala","Relli","Yadava"], "westeurasia":["Basque","English","Finnish","French", \
"Greek","Icelandic","Bergamo","Spanish","Tuscan"]}
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
        for line in f:
            spline=line.split(sep="\t")
            popName=spline[2]
            popSpline = popName.split(sep=":")
            pop = popSpline[0]
            if pop in pops[region] and lineCount%2 ==0:
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
                outCols=[str(rep), str(ss),region,archaic,str(tractLen)]
                outLine = ",".join(outCols)+"\n"
                h.write(outLine)
h.close()