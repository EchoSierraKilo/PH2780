##################################################################################
#### Part 1 LD calculation
## Q1
lab2.ped <- read.table("C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab2/lab2.ped", header=F)
dim(lab2.ped)
lab2.ped[1:5,1:10]
table(lab2.ped$V1==lab2.ped$V2)

lab2.map <- read.table("C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab2/lab2.map", header=F)
table(lab2.map$V1)

## Q3
lab2.map.cetp <- subset(lab2.map, V4 >= 56961950 & V4 <= 56983845)
cetp.snps <- lab2.map.cetp$V2
write.table(cetp.snps, file="C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab2/cetp.snps", row.names=F, col.names=F, quote=F)

## Q4 PLINK commands
# extract SNPs in cetp.snps
# plink --noweb --file lab2 --extract cetp.snps --recode --out cetp
# calculate MAF
# plink --noweb --file cetp --freq --out cetp_freq
# calculate pairwise LD (default r2 cut-off 0.2)
# plink --noweb --file cetp --r2 --out cetp_LD
# calculate pairwise LD for all SNPs
# plink --noweb --file cetp --r2 --ld-window-r2 0 --out cetp_LD_gt0
# calculate pairwise LD and output a matrix
# plink --noweb --file cetp --r2 --matrix  --out cetp_LD_matrix

## Q5 PLINK command
# type two haplotypes and save them to a file haplist
#cat > haplist
#* rs289715 rs289742
#* rs4783962 rs289715

# estimate haplotype frequencies for haplotypes in file haplist
# plink --noweb --file lab2 --hap haplist --hap-freq --out lab2_haps

#### Part 2 Quantitative traits
## Q1-2
lab2 <- read.table("C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab2/genopheno_lab2.txt", sep="\t", header=T, na.strings="X")
tapply(lab2$age, lab2$sex, mean, na.rm = T)
tapply(lab2$age, lab2$sex, sd, na.rm = T)
tapply(lab2$trait2, lab2$sex, mean, na.rm = T)
tapply(lab2$trait2, lab2$sex, sd, na.rm = T)
tapply(lab2$trait3, lab2$sex, mean, na.rm = T)
tapply(lab2$trait3, lab2$sex, sd, na.rm = T)
tapply(lab2$trait4, lab2$sex, mean, na.rm = T)
tapply(lab2$trait4, lab2$sex, sd, na.rm = T)
tapply(lab2$trait5, lab2$sex, mean, na.rm = T)
tapply(lab2$trait5, lab2$sex, sd, na.rm = T)
tapply(lab2$covar1, lab2$sex, mean, na.rm = T)
tapply(lab2$covar1, lab2$sex, sd, na.rm = T)
table(lab2$smoking, lab2$sex)
table(lab2$trait1, lab2$sex)
table(lab2$trait6, lab2$sex)
table(lab2$trait7, lab2$sex)
## Q3
par(mfrow=c(2,1))
hist(lab2$age, xlab="age", main="Histogram of Age", col="light blue")
qqnorm(lab2$age)
qqline(lab2$age, col = "red")

par(mfrow=c(2,1))
hist(lab2$trait2, xlab="trait2", main="Histogram of Trait2", col="light blue")
qqnorm(lab2$trait2)
qqline(lab2$trait2, col = "red")

par(mfrow=c(2,1))
hist(lab2$trait3, xlab="trait3", main="Histogram of Trait3", col="light blue")
qqnorm(lab2$trait3)
qqline(lab2$trait3, col = "red")

par(mfrow=c(2,1))
hist(lab2$trait4, xlab="trait4", main="Histogram of Trait4", col="light blue")
qqnorm(lab2$trait4)
qqline(lab2$trait4, col = "red")

par(mfrow=c(2,1))
hist(lab2$trait5, xlab="trait5", main="Histogram of Trait5", col="light blue")
qqnorm(lab2$trait5)
qqline(lab2$trait5, col = "red")

par(mfrow=c(2,1))
hist(lab2$covar1, xlab="covar", main="Histogram of Covariance", col="light blue")
qqnorm(lab2$covar1)
qqline(lab2$covar1, col = "red")

# repeat for trait2, covar1, log trait2, log covar1
par(mfrow=c(2,1))
hist(log(lab2$trait2), xlab="log(trait2)", main="Histogram of Log-Trait2", col="light blue")
qqnorm(log(lab2$trait2))
qqline(log(lab2$trait2), col = "red")

par(mfrow=c(2,1))
hist(log(lab2$trait4), xlab="log(trait4)", main="Histogram of Log-Trait4", col="light blue")
qqnorm(log(lab2$trait4))
qqline(log(lab2$trait4), col = "red")

par(mfrow=c(2,1))
hist(log(lab2$trait5), xlab="log(trait5)", main="Histogram of Log-Trait5", col="light blue")
qqnorm(log(lab2$trait5))
qqline(log(lab2$trait5), col = "red")

par(mfrow=c(2,1))
hist(log(lab2$covar1), xlab="log(covar1)", main="Histogram of Log-Covariance", col="light blue")
qqnorm(log(lab2$covar1))
qqline(log(lab2$covar1), col = "red")

lab2$log_trait4 <- log(lab2$trait4)
lab2$log_trait5 <- log(lab2$trait5)
lab2$log_covar1 <- log(lab2$covar1)

write.table(lab2, file="C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab2/merlin_lab2.csv", row.names=F, col.names=T, sep=",", quote=F, na="X")



