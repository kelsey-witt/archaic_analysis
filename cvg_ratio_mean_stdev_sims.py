"""
This script was written by Kelsey Witt Dillon in March 2021. The goal of this script
is to take tables that describe where simulated chromosomes have archaic ancestry,
calculate the archaic coverage of a subset of N individuals and determine the ratio 
of archaic ancestry coverage between East Asian and European simulated chromosomes.
The script outputs the mean and standard deviation of 100 replicates of each sample
size for each simulation run.

The input file is a set of replicates of a "segment_table.csv" file, a csv file containing 
simulated archaic segments and a binary score for whether simulated chromosomes from CHB or 
CEU have the archaic allele. The file is generated using the script msprime_simulation_nea_introgression.py.
The list of sample_sizes is set to N=1, 5, 10, 25, 50, 75, 100, 125, and 150, but can be 
adjusted in line 32.

Usage: python3 cvg_ratio_mean_stdev_sims.py <model>, where <model> represents the prefix of the 
segment_table.csv file that includes the introgression information (line 35). One example of this input 
is "nea_introgression_sim_1p0.02_2p0", where the first number represents the first pulse, 
the second the second pulse. The range of simNum (line 44) should be set from 1 to 1+ the 
total number of simulations to look at all of them.

This script generates a tab-delimited file labeled "<model>_mean_sd.txt". Each row represents 
a simulation model, a replicate of that simulation, and either the means or the standard 
deviations, and the columns represent the model name, the replicate name, whether the numbers 
represent means or standard deviations, and the values for each sample size, from largest to 
smallest.
"""
import csv
import random
import sys
import statistics
import math
import numpy as np

model = sys.argv[1] #formatted like nea_introgression_sim_1p0.02_2p0
sample_sizes = [1,5,10,25,50,75,100,125,150]
#sample_sizes = [1]
populations = ["CHB", "CEU"]
Pop_Tracker = {}
Pop_Cvg = {}

outfile = model + "_mean_sd.txt"
g=open(outfile, 'w')
for simNum in range(1,11):
    simFile = model + "_" + str(simNum) + "_segment_table.csv"
    cvgRatios = []
    ssMeans=[]
    ssStDev=[]
    for ss in sample_sizes:
        for repNum in range(1,101):
            for pop in populations:
                Pop_Cvg[pop] = 0
                Pop_Tracker[pop] = []
            with open(simFile) as f:
                for line in f:
                    if 'position' in line:
                        spline = line.split(sep=",")
                        for colnum in range(1,len(spline)):
                            indPop = spline[colnum][0:3]
                            Pop_Tracker[indPop].append(colnum)
                        for pop in Pop_Tracker.keys():
                            ind_downsample = random.choices(Pop_Tracker[pop],k=ss*2)
                            ind_downsample.sort()
                            Pop_Tracker[pop] = ind_downsample
                    else:
                        spline = line.split(sep=",")
                        interval = spline[0].split(sep="-")
                        intervalLen = float(interval[1])-float(interval[0])
                        for pop in populations:
                            numSites = 0
                            hasArchAllele = False
                            for i in range(0,len(Pop_Tracker[pop]),2):
                                chrom1 = spline[Pop_Tracker[pop][i]]
                                chrom2 = spline[Pop_Tracker[pop][i+1]]
                                for chrom in [chrom1,chrom2]:
                                    if int(chrom) >= 0:
                                        hasArchAllele = True
                            if hasArchAllele:
                                Pop_Cvg[pop] += intervalLen
            chb_cvg = Pop_Cvg["CHB"]
            ceu_cvg = Pop_Cvg["CEU"]
            if ceu_cvg == 0:
                cvgRatio = 0
            else:
                cvgRatio = float(chb_cvg/ceu_cvg)
            cvgRatios.append(cvgRatio)
        meanCvg = statistics.mean(cvgRatios)
        stdevCvg = statistics.stdev(cvgRatios)
        ssMeans.append(meanCvg)
        ssStDev.append(stdevCvg)
    meanLine = model + "\t" + str(simNum) + "\t mean \t" + "\t".join(str(x) for x in ssMeans) + "\n"
    g.write(meanLine)
    stDevLine = model + "\t" + str(simNum) + "\t stdev\t" + "\t".join(str(x) for x in ssStDev) + "\n"
    g.write(stDevLine)
g.close()