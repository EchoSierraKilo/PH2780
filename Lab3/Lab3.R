##################################################################################
#### Part 1 linkage of quantitative traits
## data preparation in R
lab3 <- read.table("C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab3/phenogeno_lab3.txt", sep="\t", header=T)

## trait 5 recode
lab3$trait5_new <- ifelse(is.na(lab3$trait5), "X", lab3$trait5)

## sort lab3 by famid
lab3.sorted <- lab3[order(lab3$famid),]

## save your lab3 dataset
write.table(lab3.sorted, file="C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab3/lab3.csv", row.names=F, col.names=T, sep=",", quote=F)

## work under UNIX
# check the number of columns in a file
# awk -F, '{if(NR==1) print NF}' lab3.csv
# awk -F, '{print $2,$1,$3,$4,$5,$28,$14,$15,$17,$18,$19,$20,$21,$22}' lab3.csv | sed 's/_/ /g' | sed -e '1d' > lab3.ped
# pedstats -d lab3.dat -p lab3.ped > lab3.pedstat
# merlin -d lab3.dat -p lab3.ped -m lab3.map --vc --pdf --prefix merlin_sglpt > lab3.sglpt
# merlin -d lab3.dat -p lab3.ped -m lab3.map --vc --grid 1 --pdf --prefix merlin_mltpt > lab3.mltpt
# merlin -d lab3.dat -p lab3.ped -m lab3.map --vc --grid 1 --useCovariates --pdf --prefix merlin_mltadj > lab3.mltadj

#### Part 2 Association Analysis ####
## data preparation in R
lab3 <- read.table("C:/Users/schwa/OneDrive/Documents/PH 2780 Genetic Epidemiology/Lab3/phenogeno_lab3.txt", sep="\t", header=T, na.strings=c("NA",""))
## Q1
mytable.snp1 <- table(lab3$generate, lab3$SNP1)
prop.table(mytable.snp1, margin=1) # by row
mytable.snp2 <- table(lab3$generate, lab3$SNP2)
prop.table(mytable.snp2, margin=1)
mytable.snp3 <- table(lab3$generate, lab3$SNP3)
prop.table(mytable.snp3, margin=1)
mytable.snp4 <- table(lab3$generate, lab3$SNP4)
prop.table(mytable.snp4, margin=1)
mytable.snp5 <- table(lab3$generate, lab3$SNP5)
prop.table(mytable.snp5, margin=1)

## Q2
lab3$SNP1_add <- NA
lab3$SNP1_add <- ifelse(lab3$SNP1 == "GG", 0, lab3$SNP1_add)
lab3$SNP1_add <- ifelse(lab3$SNP1 == "AG", 1, lab3$SNP1_add)
lab3$SNP1_add <- ifelse(lab3$SNP1 == "AA", 2, lab3$SNP1_add)
lab3$SNP2_add <- NA
lab3$SNP2_add <- ifelse(lab3$SNP2 == "AA", 0, lab3$SNP2_add)
lab3$SNP2_add <- ifelse(lab3$SNP2 == "AG", 1, lab3$SNP2_add)
lab3$SNP2_add <- ifelse(lab3$SNP2 == "GG", 2, lab3$SNP2_add)
lab3$SNP3_add <- NA
lab3$SNP3_add <- ifelse(lab3$SNP3 == "CC", 0, lab3$SNP3_add)
lab3$SNP3_add <- ifelse(lab3$SNP3 == "CG", 1, lab3$SNP3_add)
lab3$SNP3_add <- ifelse(lab3$SNP3 == "GG", 2, lab3$SNP3_add)
lab3$SNP4_add <- NA
lab3$SNP4_add <- ifelse(lab3$SNP4 == "GG", 0, lab3$SNP4_add)
lab3$SNP4_add <- ifelse(lab3$SNP4 == "AG", 1, lab3$SNP4_add)
lab3$SNP4_add <- ifelse(lab3$SNP4 == "AA", 2, lab3$SNP4_add)
lab3$SNP5_add <- NA
lab3$SNP5_add <- ifelse(lab3$SNP5 == "AA", 0, lab3$SNP5_add)
lab3$SNP5_add <- ifelse(lab3$SNP5 == "AG", 1, lab3$SNP5_add)
lab3$SNP5_add <- ifelse(lab3$SNP5 == "GG", 2, lab3$SNP5_add)

table(lab3$SNP1,lab3$SNP1_add)
table(lab3$SNP2,lab3$SNP2_add)
table(lab3$SNP3,lab3$SNP3_add)
table(lab3$SNP4,lab3$SNP4_add)
table(lab3$SNP5,lab3$SNP5_add)

## Q3
library(dplyr)
library(tidyverse)
# Create table for Trait 6
trait6_table <- lab3 %>%
  filter(trait6 %in% c(0, 1)) %>%
  group_by(generate, trait6) %>%
  summarise(count = n()) %>%
  spread(trait6, count) %>%
  mutate(Trait = "Trait6", 
         Control_to_Case_Ratio = `0` / `1`) 

# Create table for Trait 7
trait7_table <- lab3 %>%
  filter(trait7 %in% c(0, 1)) %>%
  group_by(generate, trait7) %>%
  summarise(count = n()) %>%
  spread(trait7, count) %>%
  mutate(Trait = "Trait7", 
         Control_to_Case_Ratio = `0` / `1`) 

# Combine tables
combined_table <- rbind(trait6_table, trait7_table)

## Q4, ANOVA
lab3$trait5_log <- log(lab3$trait5)
lab3.gp <- subset(lab3, generate=="GP")
lab3.child <- subset(lab3, generate=="CHILD")
lab3.parent <- subset(lab3, generate=="PARENT")

summary(aov(trait5_log ~ SNP3, data = lab3.gp))
bartlett.test(trait5_log ~ SNP3, data = lab3.gp)
summary(aov(trait5_log ~ SNP3, data = lab3.parent))
bartlett.test(trait5_log ~ SNP3, data = lab3.parent)
summary(aov(trait5_log ~ SNP3, data = lab3.child))
bartlett.test(trait5_log ~ SNP3, data = lab3.child)

## Q5 LM
summary(lm(trait5_log ~ SNP3, data = lab3.gp))
summary(lm(trait5_log ~ SNP3, data = lab3.parent))
summary(lm(trait5_log ~ SNP3, data = lab3.child))

## Q6-7 chi-square
#GPs
table(lab3.gp$trait7, lab3.gp$SNP4)
chisq.test(table(lab3.gp$trait7, lab3.gp$SNP4))
fisher.test(table(lab3.gp$trait7, lab3.gp$SNP4))
#Parents
table(lab3.parent$trait7, lab3.parent$SNP4)
chisq.test(table(lab3.parent$trait7, lab3.parent$SNP4))
fisher.test(table(lab3.parent$trait7, lab3.parent$SNP4))
#Children
table(lab3.child$trait7, lab3.child$SNP4)
chisq.test(table(lab3.child$trait7, lab3.child$SNP4))
fisher.test(table(lab3.child$trait7, lab3.child$SNP4))

## Q8 logistic model
summary(glm(trait7 ~ SNP3_add, family=binomial(link="logit"), data=lab3.gp))
chisq.test(table(lab3.gp$trait7, lab3.gp$SNP3_add))

## Q9 Genotypic model
# recessive model
lab3.gp$SNP3_rec <- NA
lab3.gp$SNP3_rec <- ifelse(lab3.gp$SNP3 == "CC" | lab3.gp$SNP3 == "CG", 0, lab3.gp$SNP3_rec)
lab3.gp$SNP3_rec <- ifelse(lab3.gp$SNP3 == "GG", 1, lab3.gp$SNP3_rec)
# table(lab3.gp$SNP3, lab3.gp$SNP3_rec) ## check coding
summary(glm(trait7 ~ SNP3_rec, family=binomial(link="logit"), data=lab3.gp))
chisq.test(table(lab3.gp$trait7, lab3.gp$SNP3_rec))

# dominant model
lab3.gp$SNP3_dom <- NA
lab3.gp$SNP3_dom <- ifelse(lab3.gp$SNP3 == "CC", 0, lab3.gp$SNP3_dom)
lab3.gp$SNP3_dom <- ifelse(lab3.gp$SNP3 == "CG" | lab3.gp$SNP3 == "GG", 1, lab3.gp$SNP3_dom)
# table(lab3.gp$SNP3, lab3.gp$SNP3_dom) ## check coding
summary(glm(trait7 ~ SNP3_dom, family=binomial(link="logit"), data=lab3.gp))
chisq.test(table(lab3.gp$trait7, lab3.gp$SNP3_dom))

## Q10 adjusting cov
summary(glm(trait7 ~ SNP3_add + age + sex + covar1, family=binomial(link="logit"), data=lab3.gp))
exp(0.30657)  # OR for adjusted model

summary(glm(trait7 ~ SNP3_add, family=binomial(link="logit"), data=lab3.gp))
exp(0.3714)  # OR for unadjusted model
(1.449763-1.358757)/1.358757
