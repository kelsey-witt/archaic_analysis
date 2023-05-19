"""
This script was written by Kelsey Witt Dillon in December 2020. The goal of this script
is to take a table of introgression data from simulations and determine the amount of 
archaic ancestry represented in a subset of individuals from simulated East Asian and 
European chromosomes. Sample sizes range from 1-150, and the script samples 100 replicates 
for each sample size.

The input file is a set of replicates of a "segment_table.csv" file, a csv file containing 
simulated archaic segments and a binary score for whether simulated chromosomes from CHB or 
CEU have the archaic allele. The file is generated using the script msprime_simulation_nea_introgression.py.
The list of sample_sizes is set to N=1, 5, 10, 25, 50, 75, 100, 125, and 150, but can be 
adjusted in line 31.

Usage: python3 sim_cvg_ss_multirun.py <infile>, where <infile> represents the prefix of the 
segment_table.csv file that includes the introgression information. One example of this input 
is "nea_introgression_sim_1p0.02_2p0", where the first number represents the first pulse, 
the second the second pulse. The range of simNum (line 43) should be set from 1 to 1+ the 
total number of simulations to look at all of them.

This script generates a csv file for each sample size, labeled "<infile>_[N]ind_100x.csv", 
where N is the sample size. Each row is a replicate and each column is a population, and 
the number of unique archaic sites found across the sample size of that number of individuals 
from a given population is listed in the table.
"""
import csv
import random
import sys

infile = sys.argv[1] #formatted like archaic_kfav_sim_r1_2p
sample_sizes = [1,5,10,25,50,75,100,125,150]
#sample_sizes = [1]
populations = ["CHB", "CEU"]
Pop_Tracker = {}
Pop_Cvg = {}

for ss in sample_sizes:
    outfile = infile + "_" + str(ss) + "ind_100x.csv"
    with open(outfile, 'w', newline = '') as csvfile:
        w = csv.writer(csvfile, delimiter=",")
        header = ["run"] + populations
        w.writerow(header)
        for simNum in range(1,11):
            simFile = infile + "_" + str(simNum) + "_segment_table.csv"
            for pop in populations:
                Pop_Cvg[pop] = 0
                Pop_Tracker[pop] = []
            for repNum in range(1,11):
                with open(simFile) as f:
                    for line in f:
                        if 'Chromosome' in line:
                            spline = line.split(sep=",")
                            for colnum in range(1,len(spline)):
                                indPop = spline[colnum][0:3]
                                Pop_Tracker[indPop].append(colnum)
                            for pop in Pop_Tracker.keys():
                                ind_downsample = random.sample(Pop_Tracker[pop],ss)
                                ind_downsample.sort()
                                Pop_Tracker[pop] = ind_downsample
                        else:
                            spline = line.split(sep=",")
                            for pop in populations:
                                hasArchAllele = False
                                for ind in Pop_Tracker[pop]:
                                    if "1" in spline[ind]:
                                        hasArchAllele = True
                                if hasArchAllele:
                                    Pop_Cvg[pop] += 1
                runID = str(simNum) + "-" + str(repNum)
                outline = [runID]
                for pop in populations:
                    outline.append(str(Pop_Cvg[pop]))
                w.writerow(outline)

                        

