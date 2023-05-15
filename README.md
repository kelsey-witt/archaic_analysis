# archaic_analysis

This repository contains the code needed to repeat the analyses used in "Apportioning archaic variants among modern populations", published in *Philosophical Transactions B* by Witt and colleagues in 2022. The publication is open access and available [here] (https://royalsocietypublishing.org/doi/10.1098/rstb.2020.0411).

The aim of this study was to examine how archaic ancestry is distributed across modern global populations and to look for differences in allele frequencies and the number of alleles present in each population. We compared the 1000 Genomes Phase 3 dataset to the four high-coverage archaic humans and the primary metric we used is "archaic ancestry coverage", which compares how many archaic SNPs found in an individual compares to the total number of archaic SNPs found across any individual in a population.

## Data Used
* [1000 Genomes Project Phase 3 VCF files](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/)
* High-Coverage archaic genomes for the [Altai Neanderthal](http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/), [Vindija Neanderthal](http://ftp.eva.mpg.de/neandertal/Vindija/) Chagyrskaya Neanderthal(http://ftp.eva.mpg.de/neandertal/Chagyrskaya/), and [Denisovan](http://cdna.eva.mpg.de/denisova/)
* [Simons Genome Diversity Project Data](https://reichdata.hms.harvard.edu/pub/datasets/sgdp/) - specifically the Papuan genomes
The following datasets are only needed for follow-up analyses or confirmations of our results using other datasets.
* [Sprime archaic site calls](https://data.mendeley.com/datasets/y7hyt83vxr/1)
* [Introgression tracts from Steinruecken et al. 2018](https://onlinelibrary-wiley-com.revproxy.brown.edu/doi/full/10.1111/mec.14565)
* [Introgression tracts from Sankararaman et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27032491/)

## Analyses and scripts
### Calculating allele counts, frequencies, and coverage for the 1000 Genomes Dataset
ID_rare_AFR_genotypes.py makes a csv file containing all alleles in the 1000 Genomes dataset that are low frequency (<1%) in Africa and includes the genotypes of all individuals from non-African populations.
* Usage: python3 ID_rare_AFR_genotypes.py
* Inputs: 1000 Genomes vcf and panel file. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 24), outfile (output file name, line 110), and modern_file (vcf file, lines 107 and 109)
* Output: A csv file for each chromosome where each row is a site that has an allele that is rare in African populations and the columns include position information, reference/alternative/and rare allele identities, and genotypes for all non-African individuals

ID_rare_AFR_archaic_genotypes.py takes the above output file and compares it to Neanderthal and Denisovan genomes to make a "shortlist" csv file that is limited to all alleles that are rare in Africa and shared with archaic humans. It includes the genotypes for archaic humans and non-African individuals.
* Usage: python3 ID_rare_AFR_archaic_genotypes.py
* Inputs: The csv files with genotype data for alleles that are rare in Africa (see ID_rare_AFR_genotypes.py) and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "afrInfile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in lines 62-66)
* Output: A csv file where each row is a site that has an allele that is rare in African populations and shared with archaic humans, and the columns include position information, reference/alternative/and rare allele identities, and genotypes for all archaic humans and non-African individuals

rareAFR_allele_counter.py counts archaic alleles per individual and allele frequencies of both archaic and non-archaic alleles that are rare in Africa.
* Usage: python3 rareAFR_allele_counter.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are na_only, nc_only, nv_only, d_only, n_only, n_shared, n_all, d_all, arch_shared, and arch_all, and are explained fully at the top of the script.
* Inputs: The csv files with genotype data for alleles that are rare in Africa (see ID_rare_AFR_genotypes.py) and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "afrInfile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in lines 160-164)
* Outputs: Two files: a) A csv file with counts of archaic and non-archaic alleles that are common (>20% frequency) and rare (>1% and <20%) for each population and b) A csv file with counts of archaic alleles per individual

ID_rare_AFR_snps.py makes a csv file containing all alleles in the 1000 Genomes dataset that are low frequency (<1%) in Africa and specifies whether individuals from non-African populations have that allele.
* Usage: python3 ID_rare_AFR_genotypes.py
* Inputs: 1000 Genomes vcf and panel file. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 27), outfile (output file name, line 114), and modern_file (vcf file, lines 110 and 112)
* Output: A csv file for each chromosome where each row is a site that has an allele that is rare in African populations and the columns include position information, reference/alternative/and rare allele identities, and presence/absence of the rare allele in the individuals from non-African populations (0 if absent, 1 if present)

ID_rare_AFR_archaic_snps.py takes the above output file and compares it to Neanderthal and Denisovan genomes to make a "shortlist" csv file that is limited to all alleles that are rare in Africa and shared with archaic humans. It includes the genotypes for archaic humans and specifies whether individuals from non-African populations have that allele.
* Usage: python3 ID_rare_AFR_archaic_snps.py
* Inputs: The csv files with presence/absence data for alleles that are rare in Africa (see ID_rare_AFR_snps.py) and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "afrInfile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in lines 62-66)
* Output: A csv file where each row is a site that has an allele that is rare in African populations and shared with archaic humans, and the columns include position information, reference/alternative/and rare allele identities, genotypes for all archaic humans and presence/absence of the allele in non-African individuals

archaic_mp_cvg_ss.py counts the number of unique archaic alleles found across a sample of individuals from a population and outputs 100 replicates of those counts for the sample sizes specified.
* Usage: python3 archaic_mp_cvg_ss.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are na_only, nc_only, nv_only, d_only, n_only, n_shared, n_all, d_all, arch_shared, and arch_all, and are explained fully at the top of the script.
* Input: The csv file with presence/absence data for archaic alleles that are rare in Africa (see ID_rare_AFR_archaic_snps.py). The input file path is specified in line 27 if changes are needed, and the set of sample sizes is specified in line 35.
* Output: A csv for each sample size that specifies the number of archaic alleles identified in a sample of N individuals from each population. Each sample is repeated 100 times.

Plot_genome_coverage.R is R code that is useful for taking the output of the above script and plotting it, replicating Figure 2B/2D/2F from the manuscript. As written, it will plot the coverage plot but will not export the file, so you can execute it in R studio or add commands that print the figure to a file. File names may need to be changed.

### Analyzing genome counts and coverage with the inclusion of SGDP Papuans
archaic_allele_counter_mp_papuans_n15.py counts archaic alleles per individual and allele frequencies of both archaic and non-archaic alleles that are rare in Africa.
* Usage: python3 archaic_allele_counter_mp_papuans_n15.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are d_only, n_only, and arch_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcfs and panel file, vcf files containing the SGDP Papuan genotypes per chromosome, and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "pop_file", "kgFile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in line 47 and 243-248)
* Outputs: Two files: a) A csv file with counts of archaic and non-archaic alleles that are common (>20% frequency) and rare (>1% and <20%) for each population and b) A csv file with counts of archaic alleles per individual

Allelecvg_make_table_v2_png.py makes a csv file of all alleles that are rare in Africa, present in 1000 Genomes populations and SGDP Papuans, and shared with archaic humans. It includes the genotypes for archaic humans and specifies whether individuals from non-African populations have that allele.
* Usage: python3 Allelecvg_make_table_v2_png.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are d_only, n_only, and arch_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcfs and panel file, vcf files containing the SGDP Papuan genotypes per chromosome, and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "pop_file", "kgFile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in line 41 and 194-198)
* Output: A csv file where each row is a site that has an allele that is rare in African populations and shared with archaic humans, and the columns include position information, reference/alternative/and rare allele identities, genotypes for all archaic humans and presence/absence of the allele in non-African individuals

archaic_png_cvg_ss.py counts the number of unique archaic alleles found across a sample of individuals from a population and outputs 100 replicates of those counts for the sample sizes specified.
* Usage: python3 archaic_png_cvg_ss.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are d_only, n_only, and arch_all, and are explained fully at the top of the script.
* Input: The csv file with presence/absence data for archaic alleles that are rare in Africa that includes Papuans (see Allelecvg_make_table_v2_png.py). The input file path is specified in line 31 if changes are needed, and the set of sample sizes is specified in line 32.
* Output: A csv for each sample size that specifies the number of archaic alleles identified in a sample of N individuals from each population. Each sample is repeated 100 times. The results can be plotted using Plot_genome_coverage.R (with adjustments to add PNG as a population in the legend)

archaic_venn_1kg_png.py examines how archaic variants are shared between global regions (East Asia, Europe, South Asia, and Papua New Guinea) and produces an output that can be used to generate a Venn Diagram-like figure.
* Usage: python3 archaic_venn_1kg_png.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are d_only, n_only, and arch_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcfs and panel file, vcf files containing the SGDP Papuan genotypes per chromosome, and the individual VCF files for the Altai, Chagyrskaya, and Vindija Neanderthals and the Denisovan. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are "pop_file", "kgFile", "aNeaInfile", "cNeaInfile", "vNeaInfile", and "denInfile", in line 47 and 243-248)
* Output: A csv that includes list of allele sharing patterns (a 4-digit binary number where the digits represent East Asians, Europeans, South Asians, and Papuans and 0/1 represents absence/presence in that population) and counts of the occurrences of each pattern

sank_interval_summing_papuans.py looks at archaic ancestry tract data from Sankararaman et al. 2016 and uses it to plot increases in archaic coverage with increasing sample size in global populations (see Supplementary Figure 14)
* Usage: python3 sank_interval_summing_papuans.py
* Inputs: .ind and .haplotypes files generated from Sankararaman et al. 2016. Input file paths are specified in the script itself and will need to be altered to fit your organization system. (Variables that need to be checked are indFile in line 45 and hapFile in line 69)
* Output: A csv file with the summed total of archaic tract lengths in a given sample size of individuals from "eastasia", "westeurasia", "oceania" or "southasia". The file contains all regions, all sample sizes, and both neanderthal- and denisovan-specific tracts, with 100 replicates of each population/sample/archaic source combination.

Citation: Witt KE, Villanea F, Loughran E, Zhang X and Huerta-Sanchez E. Apportioning archaic variants among modern populations. *Phil. Trans. R. Soc. B.* 377:20200411.

Please contact Kelsey Witt (kelsey_witt_dillon@brown.edu) with any questions.
