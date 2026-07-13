#R/4.5.0
dat<-read.table("minimal_pvalues.txt", as.is=T)
n_alleles = 332 ##replace with total number of tested alleles
0.05/quantile(dat$V1, probs=0.05)/n_alleles
