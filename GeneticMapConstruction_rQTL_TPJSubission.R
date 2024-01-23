#Genetic Map Construction Using R/qtl
#convert bin map in physical distances to genetic map in cM

install.packages("qtl")
library(qtl)

rils <- read.cross("csvs", genfile = "Map.csv", phefile = "Phenos.csv", estimate.map = FALSE, na.strings = "NA")
#The phefile contains the LSM data (for either aboveground [AG] or RSA data)
#mydata <- read.cross("csv", "~/Desktop", "filename.csv", na.strings="NA", genotypes=c("A","H","B"), alleles=c("A","B"), estimate.map = FALSE)
#filename = read.cross(format = "csv", "~/Desktop", "mapfile.csv", "phenofile.csv", estimate.map = F)
#estimatemap = F means that inter-marker distances will not be estimated 

#convert to selfing RILs
mydata <- convert2riself(rils)
#shows info about the data file
summary(mydata) #have 4254 markers at this point 
geno.image(mydata)
#can make a vector of colors to change how the genotypes are represented above 
#start with colors for missing genotypes, then for each of the desired genotypes (AA, BB, AB)
#cols = c("color1", etc.)

##Data Checking
#look at pattern of missing data
plotMissing(mydata)
#horizontal lines = individuals missing data (black pixels)
#vertical lines = markers missing data (black pixels)

#plot the number of genotyped markers for each individual 
#Step 1: set margins
par(mfrow=c(1,2), las=1)
#Step 2: Plot Number of Genotyped Markers 
plot(ntyped(mydata), xlab="RIL Individual", ylab = "No. Typed Markers", main = "No. Genotypes by Individual")
#Step 3: Plot Number of Genotyped Individuals 
plot(ntyped(mydata,"mar"), xlab = "Marker ID", ylab="No. Typed Individuals", main = "No. Genotypes by Marker")
#ntyped provides the numbers of genotyped markers for each individual (or the number of genotyped individuals for each marker)
#left figure = individuals missing genotypes
#right figure = markers missing genotypes

#plot MISSING data (does opposite of above but for same effect)
#missing markers
plot(nmissing(mydata, what = "mar"))
#missing individuals
plot(nmissing(mydata, what = "ind"))

#To omit the markers with lots of missing data, we first need to identify the names of the markers. We then use drop.markers().
nt.bymar <- ntyped(mydata, "mar")
todrop <- names(nt.bymar[nt.bymar < 106])
todrop
#lists marker(s) with missing data
#drop markers that are missing data in less than 106 typed individuals
#1; S06_43013565
mydata <- drop.markers(mydata, todrop) 
#have 4253 markers at this point

##Phenotype data
#change margins
par(mar=c(1,1,1,1))
#plot phenotype (choose desired phenotype)
plotPheno(mydata, pheno.col = "Height_Ctrl")
#alternatively (start with column 2 since column 1 has RIL id's)
plotPheno(rils, pheno.col = 2)

##Look for duplicate markers 
dup <- findDupMarkers(mydata, exact.only=FALSE)
#we have 1102 duplicates 
#remove duplicates
mydata <- drop.markers(mydata, unlist(dup)) 
#we are left with 2180 markers
summary(mydata)

##Check individual genotyping frequencies
g <- pull.geno(mydata)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,2), las=1)
for(i in 1:2)
  plot(gfreq[i,], ylab="Genotype Frequency", main=c("AA", "BB")[i], 
       ylim=c(0,1))
#looks good 

##Test for linkage between markers 
mydata= est.rf(mydata)
checkAlleles(mydata, threshold = 5)
#looks good (get "no apparent problems")
#Estimation of Recombination fraction
mydata = est.rf(mydata, maxit = 10000, tol = 1e-04)
#maxit = max number of iterations for EM algorithm
#EM algorithm finds max likelihood 
#tol = tolerance for determining convergence

#Identify whether markers should be dropped by investigating the chromosome length and change
#in log likelihood when one marker at a time is dropped (done per chromosome)
dropone= droponemarker(mydata, error.prob = 0.005)
#takes time, done for each chromosome
#plot results 
par(mfrow=c(2,1))
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in Chromosome Length")
#top graph shows if there is one marker that, when removed, results in an increased likelihood (i.e. above the top line, we have none)
#bottom graph shows if there are markers that give an appreciable decrease in chromosome length when omitted (each line is a decrease of 5cM)

#identify the marker on each chromosome whose omission results in the largest decreases in chromosome length as follows.
#dont drop terminal markers, check for linkage with other markers 
#drop markers 
badmar <- rownames(summary(dropone, lod.column=2))
badmar
mydata <- drop.markers(mydata, badmar) #2170 markers at this point

#re-estimate map using Kosambi
#this map function adjusts the map distance based on interference which 
#changes the proportion of double crossovers and converts genetic distances into rf's
newmap= est.map(mydata, map.function = "kosambi", n.cluster = 10, tol = 1e-04, maxit  = 1000)
#rename the map here to keep original file as class "cross"
#n.cluster = number of parallel nodes to use 
#takes a moment

par(mfrow=c(1,1), las=1)
plotMap(newmap, show.marker.names = F)
#do not show marker names on map; makes it hard to visualize 

#save Genetic Map as .tiff file
plotMap(mydata2)
tiff(filename = "mydata_GeneticMap.tiff", width = 11, height = 8.5, units = "in", res = 900)
plotMap(mydata2)
dev.off()

#replace old map with new map (if desired)
mydata=replace.map(mydata, newmap)
summaryMap(mydata)

##Look for problem individuals (not a crucial step)
#investigate observed number of crossovers in each individual
#count number of crossovers
xo = countXO(ril_ht2)
plot(xo, ylab = "No. Crossovers")
#shows number of individuals with a certain number of crossovers (many of which have more than 20)
xo[xo>20]
#get list of RILs with more than 20 XO's (gives # of XO's)

#identifying double-crossovers by calculating genotyping error lod scores 
#Estimate genotyping error rate from the data as a function of est.map
#this estimates the inter-marker distances WHILE ALSO calculating the log likelihood for each chromosome
#If we run est.map with different values for error rate (i.e. error.prob) we can identify the
#maximum likelihood estimate of the error rate
loglik=err= c(0.001, 0.0025, 0.005, 0.0075, 0.010, 0.0125, 0.015, 0.0175, 0.02)
for (i in seq(along=err)) {
  cat(i, "of", length(err), "\n")
  tempmap= est.map(map.rils, error.prob = err[i])
  loglik[i] = sum(sapply(tempmap, attr, "loglik"))}

lod= (loglik - max(loglik))/log(10)

plot(err, lod, xlab="Genotyping Error Rate", xlim = c(0, 0.02), ylab = expression(paste(log[10], "likelihood")))
#want close to 0
#the lod score compares the likelihood for a genotype being an error versus not being an error 
#ours is very close to 0

#Identify Double Crossovers 
rilmap2 <- calc.errorlod(rilmap2, error.prob=0.005) #takes some time
toperr= top.errorlod(rilmap2, cutoff = 1.45)
#cutoff = 5 in manual,  will be more stringent here since we have a larger 
#data set and lower sequencing depth on the RILs
print(toperr)
write.csv(toperr, "~/Desktop/TopErr.csv", quote = F, row.names = F)
#export toperr file as csv onto desktop; makes it easier to compare to genetic map
#look at plot of flagged genotypes on a particular chromosome
plotGeno(map.rils, chr=1, ind=toperr$id[toperr$chr==1],cutoff = 3, include.xo = F)
#flagged cases are when markers switch from A to B to A (or vice versa)
#exclude H's here since hets they were excluded in the conversion to RIL cross

#delete specific markers that appear to be double crossovers upon investigation 
mydata$geno[[7]]$data[mydata$pheno$id=="N65_F5","S07_14576705"]= NA
mydata$geno[[7]]$data[mydata$pheno$id=="E18-W_F5","S07_14576705"]=NA
mydata$geno[[6]]$data[mydata$pheno$id=="E30_F5","S06_12833783"]=NA
mydata$geno[[1]]$data[mydata$pheno$id=="E40-W_F5","S01_433935"]=NA
mydata$geno[[5]]$data[mydata$pheno$id=="E21-W_F5","S05_12497285"]=NA
mydata$geno[[9]]$data[mydata$pheno$id=="N101_F4","S09_23571694"]=NA
mydata$geno[[7]]$data[mydata$pheno$id=="N111_F5","S07_13013790"]=NA
mydata$geno[[7]]$data[mydata$pheno$id=="N182_F5","S07_59004230"]=NA
mydata$geno[[5]]$data[mydata$pheno$id=="L17_F5","S05_78366"]=NA
mydata$geno[[7]]$data[mydata$pheno$id=="E15_F5","S07_65410809"]=NA
mydata$geno[[5]]$data[mydata$pheno$id=="E8_F5","S05_71739949"]=NA
mydata$geno[[4]]$data[mydata$pheno$id=="N124_F5","S04_68562406"]=NA
mydata$geno[[10]]$data[mydata$pheno$id=="L37_F5","S10_67696"]=NA

summary(map.rils)

#Plot final map 
mydata.clean= mydata
plotMap(mydata.clean)
tiff(filename = "GeneticMap.tiff", width = 11, height = 8.5, units = "in", res = 900)
plotMap(mydata.clean)
dev.off()
summaryMap(mydata.clean)
summary(mydata.clean)

