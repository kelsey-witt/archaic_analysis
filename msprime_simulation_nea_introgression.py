"""
This script was written by Fernando Villanea and Kelsey Witt in 2020. It is a demographic
model of Africans, East Asians and Europeans, with a Neanderthal lineage that admixes
with modern humans once or twice depending on user inputs. The model is based on the
Gravel model, modified with parameters from Moorjani et al. 2016.

Usage: python3 msprime_simulation_nea_introgression.py <firstPulse> <secondPulse>, where
firstPulse refers to the percent of admixture from Neanderthals into the ancestor of
modern Europeans and East Asians, while secondPulse refers to the percent of admixture 
from Neanderthals into the ancestor of East Asians only. Both numbers are percentages 
expressed as decimals (2% = 0.02). For a single pulse, set secondPulse to 0. This program 
requires numpy, msprime version 0.7.2, and gsl 2.5 to run.

Any of the variables can be changed depending on how you want to adjust the model. Broadly,
variables starting with N refer to population sizes, variables starting with T refer to
time-based events (ie population splits/expansions/admixture events), and variables
starting with m refer to migration rates. num_rep in line 86 sets the number of replicates 
and length of chromosome to be simulated is in line 153 in the "sims" variable. The number
of modern simulated genomes are num_<pop> in line 80-82.

This script outputs a .csv file for each replicate where the columns are each simulated
european and east asian chromosome, the rows are the segments with modeled neanderthal 
ancestry, and the individuals are scored for presence/absence of the archaic segment
for each segment observed (0 if absent, 1 if present)
"""
import math
import msprime
import numpy as np
import random
import os
import sys

#Demographic model based on Gravel model, modified with  parameters from Moorjani et al 2016
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
	N_BE = 1000
   
    # Times are provided in years, so we convert into generations.
	generation_time = 25
	T_archaic = 300e3 / generation_time
	T_AF = 220e3 / generation_time
	T_B = 140e3 / generation_time
	T_BE = 50e3 / generation_time
	T_pulse_1_start = 44.3e3/ generation_time
	T_pulse_1_end = T_pulse_1_start-20
	T_EU_AS = 40e3 / generation_time
	T_pulse_2_start = 39e3 / generation_time
	T_pulse_2_end = T_pulse_2_start-20
	T_dilution = 25e3 / generation_time
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
	m_archaic_B = float(sys.argv[1]) #0.02 <- this is where you adjust the first pulse of introgression into both
	m_archaic_AS = float(sys.argv[2]) #0.002 <- this is where you adjust the first pulse of introgression into Asia only	
	# Sample sizes for modern pops
	num_as = 206 #394
	num_eur = 198 #170
	num_afr = 204
	# Population IDs correspond to their indexes in the population
	# configuration array. Therefore, we have 3=archaic 0=YRI, 1=CEU and 2=CHB 4=dilution pop
	# initially. Archaic single sample must be last in the sampled populations for simplicity
	num_rep= 100
	mu= 1.5e-8
	rho= 1.0e-8
	samples = [msprime.Sample(population=0,time=0)]*num_afr #sample africans
	samples.extend([msprime.Sample(population=1,time=0)]*num_eur) #sample europeans
	samples.extend([msprime.Sample(population=2,time=0)]*num_as) #sample asians
	samples.extend([msprime.Sample(population=3,time=T_EU_AS)]*(1)) #sample 1 archaic for comparison
	population_configurations = [
		msprime.PopulationConfiguration(initial_size=N_AF),
		msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU),
		msprime.PopulationConfiguration(initial_size=N_AS, growth_rate=r_AS),
		msprime.PopulationConfiguration(initial_size=N_archaic),
		msprime.PopulationConfiguration(initial_size=N_BE)
	]
	migration_matrix = [
		[	0, m_AF_EU, m_AF_AS,	0,	0],
		[m_AF_EU,	0, m_EU_AS,	0,	0],
		[ m_AF_AS, m_EU_AS,	0,	0,	0],
		[	0,	0,	0,	0,	0],
		[	0,	0,	0,	0,	0],
	]
	demographic_events = [
	#Second pulse to AS
		msprime.MassMigration(time=T_pulse_2_start, source=2, destination=3, proportion=m_archaic_AS),
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
    # Basal Eurasian merge into B
		msprime.MassMigration(
                        time=T_BE, source=4, destination=1, proportion=1.0),
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

	sims = msprime.simulate(samples=samples,population_configurations=population_configurations,migration_matrix=migration_matrix,demographic_events=demographic_events,record_migrations=True,length=10e7,mutation_rate=mu,recombination_rate=rho,num_replicates=num_rep)
	print('done simulating')
	#record output 
	simNum = 101
	for sim in sims:
		print(str(simNum))
		lenOutfile = "nea_introgression_sim_1p" + str(m_archaic_B) + "_2p" + str(m_archaic_AS) + "_" + str(simNum) + "_segment_table.csv"
		g = open(lenOutfile, 'w')
		headercols = ["position"] 
		for i in range(0,num_eur):
			ind = "CEU" + str(i)
			headercols.append(ind)
		for j in range(0,num_as):
			ind = "CHB" + str(j)
			headercols.append(ind)
		headerline = ",".join(headercols)+"\n"
		g.write(headerline)
		Europeans = sim.get_samples(1) #pulls samples of a given pop ID
		EastAsians = sim.get_samples(2)
		for tree in sim.trees():
			position = tree.get_interval()
			interval = "-".join(map(str,tree.interval))
			archaicLn = [0]*(num_as+num_eur)
			cur_node = len(samples)-1  #position of the very last leaf, when adding more modern pops make sure archaic is still last
			while tree.get_time(tree.get_parent(cur_node)) < T_archaic: #this updates the neanderthal parent node to the oldest one before human-archaic split
				cur_node = tree.get_parent(cur_node) #get parent node of introgressed leaf multiple times if not the oldest one before human-archaic split
			if len(list(tree.leaves(cur_node)))>1: #this is just to speed things up, if there are only archaic leaves, it doesn't do the next part
				lenLine = [interval]
				segLength = tree.num_sites
				archaicSeg = False
				for l in tree.leaves(cur_node): #for each leaf in archaic node
					if l in Europeans or l in EastAsians: #if the leaf is in the right population
						archaicLn[(l-num_afr)]=1 #set it to 1
						if segLength>=0:
							archaicSeg = True
					elif tree.get_population(l)== 3: continue
					elif tree.get_population(l)== 4: print('oops,Dilution') #sanity check
					else:
						print('oops,Africans') #sanity check
				if archaicSeg:
					for l in archaicLn:
						lenLine.append(l)
					dataline = ",".join(map(str,lenLine))+"\n"
					g.write(dataline)
		simNum+=1
model = out_of_africa()
