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
* [Introgression tracts from Skov et al. 

## Analyses and scripts
### Calculating archaic allele density
archaic_snp_density_sprime_admixed.py divides the genomes of the admixed American populations (currently set to PEL, MXL, CLM, and PUR) into ancestry segments, and then calculates the % of sites within those ancestry segments that are classified as archaic by Sprime.
* Usage: python3 archaic_snp_density_sprime_admixed.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are nd_either, nd_both, n_only, n_all, d_only, d_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcf and panel file, Sprime calls, and modern ancestry tract calls. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 26), outfile (output file name, line 25), infile (the ancestry tract bed files, line 83), modern_file (vcf file, line 137), and sprime_file (sprime calls, line 140)
* Output: A csv file where each row is a different ancestry tract in a different individual and the columns are, in order, the population and ancestry call separated by an underscore, the individual ID, the coordinates of the tract in standard BED format, the number of archaic alleles, the length of the tract, and the density (the number of archaic alleles divided by the tract length)

Citation: Witt KE, Villanea F, Loughran E, Zhang X and Huerta-Sanchez E. Apportioning archaic variants among modern populations. *Phil. Trans. R. Soc. B.* 377:20200411.

Please contact Kelsey Witt (kelsey_witt_dillon@brown.edu) with any questions.
