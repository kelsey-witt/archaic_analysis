#This script takes data from Steinruecken et al. 2018 and generates
#a plot to show how the ratio of East Asian to European archaic
#coverage changes with increasing sample size. See Supplemental
#Figure 3 from Witt et al. 2022.

AS <- read.table("./CHBS_lax_data_2_no_denisova.bed",sep="\t",header=FALSE)
EU <- read.table("./CEU_lax_data_2_no_denisova.bed",sep="\t",header=FALSE)

#exclude chromosome and position columns
myvars <- names(AS) %in% c("V1", "V2", "V3")
AS10 <- AS[!myvars]

myvars <- names(EU) %in% c("V1", "V2", "V3") 
EU10 <- EU[!myvars]

#clean the memory
rm(AS)
rm(EU)

#Assign identity bt prob cut-off
cut_off <- 0.45

AS10[AS10<=cut_off]<-0
AS10[AS10>cut_off]<-1

EU10[EU10<=cut_off]<-0
EU10[EU10>cut_off]<-1

(mean(rowSums(AS10)/394))/(mean(rowSums(EU10)/170))
(mean(colSums(AS10)/26959))/(mean(colSums(EU10)/26959))
AS_ind_freq = colSums(AS10)/26959
EU_ind_freq = colSums(EU10)/26959

EU_ind_counts <- rowSums(EU10)
EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)

AS_sample = sample(1:394, 170, replace = F)
AS_10_sampled = AS10[,AS_sample]
AS_ind_counts <- rowSums(AS_10_sampled)
AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
AS_coverage/EU_coverage

#loop goes here for 1,5,10,25,50,75,100,125,150,170
#sampling
coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 1, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 1, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_1 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 5, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 5, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_5 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 10, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 10, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_10 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 25, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 25, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_25 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 50, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 50, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_50 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 75, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 75, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_75 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 100, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 100, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_100 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 125, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 125, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_125 = coverage_n

coverage_n = vector()
for(i in 1:100){
  EU_sample = sample(1:170, 150, replace = F)
  EU_10_sampled = EU10[,EU_sample]
  #EU_ind_counts <- rowSums(EU_10_sampled)
  EU_ind_counts <- EU_10_sampled
  EU_coverage = length(which(EU_ind_counts != 0))/length(EU_ind_counts)
  AS_sample = sample(1:394, 150, replace = F)
  AS_10_sampled = AS10[,AS_sample]
  #AS_ind_counts <- rowSums(AS_10_sampled)
  AS_ind_counts <- AS_10_sampled
  AS_coverage = length(which(AS_ind_counts != 0))/length(AS_ind_counts)
  coverage_n = c(coverage_n,(AS_coverage/EU_coverage))
}
emp_coverage_150 = coverage_n

stein_data <- data.frame(emp_coverage_1,emp_coverage_5,emp_coverage_10,emp_coverage_25,emp_coverage_50,emp_coverage_75,emp_coverage_100,emp_coverage_125,emp_coverage_150)
str(emp_coverage)
boxplot(emp_coverage, names = sample_size, xlab = "Individuals sampled (100 Replicates)", ylab = "Ratio of Genome Coverage East Asia/Europe")
abline(h=1.0, col = "red",lty=2)
abline(h=1.20028, col = "blue",lty=2)