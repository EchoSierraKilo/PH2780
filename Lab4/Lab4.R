######################################################################################
#### Part 1 GWAS ####
## Q1
map <- read.table("lab4.map", header=FALSE)
geno <- read.table("lab4.ped", header=FALSE)
table(map$V1)
dim(map)
dim(geno)

## Q2-Q10 plink and UNIX commands, see Procedure handout

## Q11 Manhattan plot, QQ plot and lambda calculation
hdl.results <- read.table("hdl.assoc.linear", header = TRUE)
hdl.results.add <- subset(hdl.results, TEST == "ADD")

## Manhattan plot
pdf("MH.pdf")
plot(hdl.results.add$BP, -log10(hdl.results.add$P), col="blue")
dev.off()

## QQ plot
qqpval = function(x){
  x <- sort(-log10(x[x>0]))
  n <- length(x)
  pp <- ppoints(n)
  plot(-log10(rev(pp)), x, xlab="Expected", ylab="Observed")
  abline(0,1,lty=2, col="purple")
}

pdf("QQ.pdf")
qqpval(hdl.results.add$P)
dev.off()

## lambda calculation
chi2 <- qchisq(hdl.results.add$P, df=1, lower.tail=F) 	
g <- median(chi2, na.rm=TRUE)		
lambda <- round(g/qchisq(0.5, df=1), digits=5)

## Q12 prepare meta-analysis data ##
hdl.results <- read.table("hdl.assoc.linear", header = TRUE)
hdl.results.add <- subset(hdl.results, TEST == "ADD")
hdl.results.add <- hdl.results.add[,c("SNP", "BP", "NMISS","BETA", "SE", "P")]
lab4.map <- read.table("lab4_frq.frq", header=TRUE)
hdl.results.add <- merge(lab4.map, hdl.results.add, by="SNP", all.x=F, all.y=T)
write.table(hdl.results.add, file="hdl_metal_file1.txt", row.names=F, col.names=T, sep=" ", quote=F)

## follow procedure to generate metal file ##

## Q12 check meta-analysis results file ##
metal.results <- read.table("hdl_metal_results1.txt", header=T)
metal.results.sort <- metal.results[order(metal.results$P.value),]
head(metal.results.sort)

## Q13 meta-analysis MH plot and QQ plot, lambda calculation

lab4.map <- read.table("lab4.map", header=FALSE)
names(lab4.map) <- c("CHR", "SNP", "cM", "BP")
metal.results <- merge(lab4.map, metal.results, by.x="SNP", by.y="MarkerName", all.x=F, all.y=T)

pdf("MH_meta.pdf")
plot(metal.results$BP, -log10(metal.results$P.value), col="blue")
dev.off()

qqpval = function(x){
  x <- sort(-log10(x[x>0]))
  n <- length(x)
  pp <- ppoints(n)
  plot(-log10(rev(pp)), x, xlab="Expected", ylab="Observed")
  abline(0,1,lty=2, col="purple")
}

pdf("QQ_meta.pdf")
qqpval(metal.results$P.value)
dev.off()

## lambda calculation
chi2 <- qchisq(metal.results$P.value, df=1, lower.tail=F) 	
g <- median(chi2, na.rm=TRUE)		
lambda <- round(g/qchisq(0.5, df=1), digits=5)

