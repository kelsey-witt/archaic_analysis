"""
This script was written by Kelsey Witt in November 2021. The goal of this script is
to generate a vcf file of simulated chromosomes as well as a list of "true"
introgressed sites to assess the performance of the SNP-counting method. The
model includes a single introgression event from Neanderthals into the ancestor of 
Europeans and East Asians.

Usage: python3 lewontin_sim_truesites.py. This program requires msprime and numpy to
run.

This script generates two output files for each replicate simulation (established in 
num_rep, line 73). The first file (nea_sim_trueSiteVSvcf_#.vcf, line 140) is a .vcf file 
containing SNP information for all simulated individuals, and is used for identifying archaic
SNPs. The second (nea_sim_trueSiteVSvcf_trueSiteList#.txt, line 141) is a .txt file that 
includes all sites that are introgressed from archaic humans and shared with the simulated 
Neanderthal individual, and is used as a list of "true sites" to compare to.
"""
import math
import msprime
import numpy as np
import random
import os

#Demographi model based on Gravel model, modified with  parameters from Moorjani et al 2016
#Ne0 Neanderthal Ne 2500
#Ne1 Eur Ne 10000
#Ne2 Asn Ne 10000
#mu 1.5e-8 bp/gen
#rho 1.0e-8 bp/gen
#Shared pulse of introgression 44,300
#Eurasian split at 40,000 based on Mellars P 2006
#Asian pulse at 39,000 or Neanderthals would be extinct

def out_of_africa():
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000 #size of bottleneck
    N_AS0 = 510 #size of bottleneck
    N_archaic = 2500
   
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_archaic = 300e3 / generation_time
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_pulse_1_start = 44.3e3/ generation_time
    T_pulse_1_end = T_pulse_1_start-20
    T_EU_AS = 40e3 / generation_time
    T_pulse_2_start = 39e3 / generation_time
    T_pulse_2_end = T_pulse_2_start-20
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004 #recovery rate (growth) after bottleneck
    r_AS = 0.0055 #recovery rate (growth) after bottleneck
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    m_archaic_B = 0.05 #0.02 <- this is where you adjust the first pulse of introgression into both
    # Sample sizes for modern pops
    num_as = 206 #394
    num_eur = 198 #170
    num_afr = 204
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 3=archaic 0=YRI, 1=CEU and 2=CHB 4=dilution pop
    # initially. Archaic single sample must be last in the sampled populations for simplicity
    num_rep= 20
    mu= 1.5e-8
    rho= 1.0e-8
    samples = [msprime.Sample(population=0,time=0)]*num_afr #sample africans
    samples.extend([msprime.Sample(population=1,time=0)]*num_eur) #sample europeans
    samples.extend([msprime.Sample(population=2,time=0)]*num_as) #sample asians
    samples.extend([msprime.Sample(population=3,time=T_EU_AS)]*(2)) #sample 2 archaic for comparison
    population_configurations = [
        msprime.PopulationConfiguration(initial_size=N_AF),
        msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(initial_size=N_AS, growth_rate=r_AS),
        msprime.PopulationConfiguration(initial_size=N_archaic)
    ]
    migration_matrix = [
        [	0, m_AF_EU, m_AF_AS,	0],
        [m_AF_EU,	0, m_EU_AS,	0],
        [ m_AF_AS, m_EU_AS,	0,	0],
        [	0,	0,	0,	0],
    ]
    demographic_events = [
    # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(0,1)),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(1,0)),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(0,2)),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(2,0)),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(1,2)),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0, matrix_index=(2,1)),
        msprime.MassMigration(
                        time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(
                        time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
    # First pulse to B
        msprime.MassMigration(time=T_pulse_1_start, source=1, destination=3, proportion=m_archaic_B),
    # Population B merges into YRI at T_B
        msprime.MigrationRateChange(time=T_B, rate=0, matrix_index=(0,1)),
        msprime.MigrationRateChange(time=T_B, rate=0, matrix_index=(1,0)),
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
    # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0),
    # Population YRI merges into A at T_archaic_AF
        msprime.MassMigration(
            time=T_archaic, source=0, destination=3, proportion=1.0),
    # Size changes to N_A
        msprime.PopulationParametersChange(
            time=T_archaic, initial_size=N_A, population_id=3)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
#	dd = msprime.DemographyDebugger(
#		population_configurations=population_configurations,
#		migration_matrix=migration_matrix,
#		demographic_events=demographic_events)
#	dd.print_history()

    sims = msprime.simulate(samples=samples,population_configurations=population_configurations,migration_matrix=migration_matrix,demographic_events=demographic_events,record_migrations=True,length=10e6,mutation_rate=mu,recombination_rate=rho,num_replicates=num_rep)
    print('done simulating')
    #record output 
    simNum = 1
    for sim in sims:
        print(str(simNum))
        #Write simulated data as a VCF of diploid individuals (the number of sampled chromosomes must ne divisible by two)
        vcfOutfile = "nea_sim_trueSiteVSvcf_" + str(simNum) + ".vcf"
        sitesOutfile = "nea_sim_trueSiteVSvcf_trueSiteList" + str(simNum) + ".txt"
        Europeans = sim.get_samples(1)
        EastAsians = sim.get_samples(2)
        with open(vcfOutfile, "w") as vcf_file:
            sim.write_vcf(vcf_file, ploidy=1)
        posList = []
        with open(vcfOutfile,'r') as h:
            for line in h:
                if '#' not in line:
                    spline = line.split()
                    position = spline[1]
                    posList.append(position)
        with open(sitesOutfile, 'w') as f:
            headCols = ["position","leaves","population","shared_with_sampled_nea"]
            headLine = "\t".join(headCols) + "\n"
            f.write(headLine)
            for mr in sim.migrations(): #for all introgression events
                if (mr.source == 1 and mr.dest == 3): # if the segment is from archaic introgression
                    for tree in sim.trees(leaf_lists=True): #for each segment of the genome, include list of samples
                        if mr.left > tree.get_interval()[0]: #the segment is before the archaic segment, keep going
                            continue
                        if mr.right <= tree.get_interval()[0]: #the segment is after the archaic segment, stop looking
                            break
                        leafList = []
                        left = tree.get_interval()[0]
                        right = tree.get_interval()[1]
                        leafList.append(mr.node)
                        cur_node = mr.node
                        inSampledNea = False
                        while tree.get_time(tree.get_parent(cur_node)) < T_archaic: #this updates the neanderthal parent node to the oldest one before human-archaic split
                            cur_node = tree.get_parent(cur_node)
                            leafList.append(cur_node)
                        inEuro = False
                        inAsia = False
                        for mutation in tree.mutations():
                            leafOut = ""
                            popStatus = ""
                            neaStatus = ""
                            if mutation.node in leafList:
                                for l in tree.leaves(mutation.node):
                                    if l in Europeans:
                                        inEuro = True
                                        leafOut += str(l) + ","
                                    elif l in EastAsians:
                                        inAsia = True
                                        leafOut += str(l) + ","
                                    elif str(l) == "608" or str(l) == "609": #if Neanderthal is in the introgressed node
                                        inSampledNea = True
                                if inEuro and inAsia:
                                    popStatus = "both"
                                elif inEuro:
                                    popStatus = "europeans"
                                elif inAsia:
                                    popStatus = "eastasians"
                                leafOut = leafOut[:-1]
                                neaStatus = "yes" if inSampledNea else "no"
                                outCols = [str(posList[mutation.site]),leafOut,popStatus,neaStatus]
                                outLine = "\t".join(outCols) + "\n"
                                f.write(outLine)
        simNum+=1
model = out_of_africa()
