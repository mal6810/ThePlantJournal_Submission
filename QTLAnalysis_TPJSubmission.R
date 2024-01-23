# QTL mapping, single QTL scans 
######
rils_84_impXO <- calc.genoprob(rils_84_impXO, step = 1, error.prob = 0.0001, map.function = "kosambi")

rils_84_impXO.sc1 <- scanone(rils_84_impXO, pheno.col = 1:9, method = "em", model = "normal")

rils_84_impXO.perms <- scanone(rils_84_impXO, pheno.col = 1:9, method = "em", model = "normal", n.perm =1000, verbose=TRUE, n.cluster = 10)

##############################################################################################
#(1) Produce a text file indicating the most significant marker on each chromosme is extracted 
#(they may or may not cross the significance threshold determined by the permutation test.), 
#the drop 1.0 LOD intervals, LODs, and pvalues for QTL detected.

##########################################################################################
sum <- summary(rils_84_impXO.sc1, threshold = 0, perms = rils_84_impXO.perms, pvalues = TRUE, format = "tabByCol", ci.function = "lodint", drop=1.0, expandtomarkers=TRUE)

space <- " "
for(i in 1:9){
  x<-capture.output(sum[[i]])
  cat(colnames(rils_84_impXO$pheno[i]), file="10_Initial QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(x, file="10_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
  cat(space, file="10_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
}

#need to omit pheno col 1 [id] and keep columns 2-11
sum_AG <- summary(AG.sc1, threshold = 0, perms = AG.perms, pvalues = TRUE, format = "tabByCol", ci.function = "lodint", drop=1.0, expandtomarkers=TRUE)
space <- " "

i = numeric(10) #do this first to make a vector. . .
for(i in 1:10){
  x_AG<-capture.output(sum_AG[[i]]) 
  cat(colnames(AG_rils$pheno[2:11]), file="10_AG_Initial QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(x_AG, file="10_AG_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
  cat(space, file="10_AG_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
}
##########################################################################################

#(2) Export csv file with LOD scores at every marker. 

##########################################################################################
lods<-as.data.frame(rils_84_impXO.sc1)
write.csv(lods, file = "11_LOD scores for every marker.csv", row.names = FALSE)

##########################################################################################
#(3) Produce PDF file indicating genome scans with alpha=.05 threshold marked
##########################################################################################
z<-summary(rils_84_impXO.perms, alpha=.05) ### this object gives threshold summary
pdf(file="12_QTL Plots1.pdf", width=11, height=8.5)
for(i in 1:9){
  plot(rils_84_impXO.sc1, lodcolumn=i, lwd=2.0, gap=0, bandcol="gray70", incl.markers=TRUE, main=colnames(rils_84_impXO$pheno[i]), xlab=c("Threshold for alpha=.05 using 1000 permutations", z[i]))
  add.threshold(rils_84_impXO.sc1, perms=rils_84_impXO.perms, alpha=0.05, lodcolumn=i, gap=0)
}
dev.off()

#for above, sc1 file only uses phenotype columns 2-11 from AG_rils file; reflect this in main=colnames(AG_rils$pheno[2:11])
#can instead/also remove "main=colnames" entirely 
##########################################################################################

#(4) Produce a csv file with the map positions for each marker.

#*NOTE*# This may produce a warning consistent with the number of chromosomes. 
#This is expected as the program is making the user aware that it is appending the column names to the file. 

##########################################################################################
newmap<-pull.map(rils_84_impXO)
for(i in 1:length(names(newmap))){
  snps<-names(newmap[[i]])
  gm<-c(snps, newmap[[i]])
  gm2<-matrix(gm, ncol=2)
  write.table(file="13_Genetic Map Positions.csv", sep=",", append=TRUE, gm2)
}

##########################################################################################

##########################################################################################

#(5)run sim.geno (necessary for some of the following functions)

##########################################################################################
cross2<-sim.geno(rils_84_impXO)
##########################################################################################

######################################################################################
########																	##########
######## Steps #(23) - #(32) represent a loop that will scan for additional ##########
######## qtl, run model selection, and fit the qtl model for each phenotype ##########
######## in turn automatically. Run the entire code (Steps #(23) - #(32) )  ##########
######## simultaneously.                                                    ##########
########																	##########
######################################################################################


######################################
##Model Selection and QTL refinement##
######################################

#(6) Specifies the phenotype you wish to work with. If you are troubleshooting or otherwise want to 
#run each step individually for a given phenotype rather than the whole things a loop, just replace 
#k in ' phenotype<-k ' with the number of the phenotype you want to use, ignore the line ' for(k in 2:nphe(cross)){ '  
#and run each step individually. 

##########################################################################################
phenotype = numeric()

for(k in 2:nphe(AG_rils)){
  phenotype = c(phenotype, k)
  pheno_AG<-colnames(AG_rils$pheno[phenotype])
  print(pheno_AG)}

##########################################################################################

#(24) Makes a qtl object containing the ALL the qtl from file="10_Initial QTL hits by phenotype" 
#for the phenotype you specified.  

##########################################################################################
#in sum_AG file, make each phenotype its own data frame
sum_AG_htc = as.data.frame(sum_AG$Height_Ctrl) #repeat for each phenotype
#then rename the first column as "Marker" (repeat for each phenotype)
sum_AG_htds <- sum_AG_htds %>% 
  rownames_to_column(var = "marker")
#add a column to the dataframe for the phenotype 
sum_AG_htc$pheno <- "Height_Ctrl" #make sure these names match the names in sum_AG
#then merge the dataframes in batches
sum_AG_all = merge(x = sum_AG_htc, y = sum_AG_htds, by = c("marker", "chr", "pos", "ci.low", "ci.high", "lod", "pval", "pheno"), all =  T, sort = F)
#do this in steps; keep adding one phenotype add a time --> final dataframe = sum_AG_allpheno

chromo<-sum[[pheno]] #this is just a list
#got the above to work by: 
chromo_AG=sum_AG[pheno_AG]
#convert chromo to dataframe since I get an error for nrows below otherwise 
chromo_AG2 = as.data.frame(chromo_AG) 
#the markers column is wrong here; need to fix (markers are the same for each phenotype, but should be unique)

chr<-{}
pos<-{}
for(i in 1:nrow(chromo)){
  chr1<-chromo[i,1]
  chr2<-as.numeric(as.character(chr1))
  chr<-c(chr, chr2)
}
for(i in 1:nrow(chromo)){
  pos1<-chromo[i,2]
  pos<-c(pos, pos1)
}
#in below, only changed rils_84_impXO to my AG_rils "cross" file
if(is.na(chr[1])==FALSE) {
  qtl<-makeqtl(rils_84_impXO, chr=chr, pos=pos, what="prob")
  createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
  formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))} 
##########################################################################################

#(25) Scans for additional linked QTL conditioning on the QTL already detected. .

#*NOTE*# You may receive warning messages about dropping individuals with missing phenotype data 
#and/or that the column names in scanone input do not match those in perms input. These are both 
#expected under some circumstances and do not effect the output of the code. Also you may receive 
#a warning that there is no chromosome number NA if no additional QTL are to be found.

##########################################################################################
rils_84_impXO.aq<-addqtl(rils_84_impXO, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", model="normal")
#for pheno.col, use the number of the phenotype (i.e. 10) or columns name and do one at a time 
#do the above for each pheno; will get warning when dropping individuals with missing phenotypes
Ht_Ctrl.aq<-addqtl(AG_rils, pheno.col="Height_Ctrl", qtl=qtl, formula=formula, method="hk", model="normal")
#repeat for each phenotype

colnames(ht_ctrl.aq)[colnames(ht_ctrl.aq) == "lod"] <- "Ht_Ctrl"
#rename the last column in each pheno_treatment.aq file so it's reflective of each phenotype
#make sure names here are the SAME as in the AG.perms file

AG_combined.aq = c(file1.aq, file2.aq, etc) #essentially a cbind function
#then combine all of the pheno_trtmnt.aq files 

#ignore 
AG_combined3_aq = AG_combined2.aq %>% select(-Marker, -chr, -pos) #this removes desired columns

sub.perms <-subset(rils_84_impXO.perms, lodcolumn=phenotype)
AG_sub.perms = subset(AG.perms, lodcolumn=pheno_AG)
#subsets LOD scores for just the phenotypes 

xx<-capture.output(summary(rils_84_impXO.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.0, expandtomarkers=TRUE))
AG_xx2<-capture.output(summary(AG_combined.aq, perms=AG_sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.0, expandtomarkers=TRUE))

xxy<-c(pheno,xx)
cat(xxy, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
cat(space, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
##########################################################################################

#(26) Adds additional QTL (if any) in the file "14_Additional QTL hits by phenotype.txt" 
#to the qtl object and update the formula object.

##########################################################################################
AG_sum.aq<-summary(AG_combined.aq, perms=AG_sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=0.95, expandtomarkers=TRUE)
#the above provides output summary of additional qtl
AG_sum.aq_agbc <- AG_sum.aq_agbc %>% 
  rownames_to_column(var = "marker")
AG_sum.aq_agbc$pheno <- "AGB_Ctrl"
AG_sum.aq_AGB = merge(x = AG_sum.aq_agbc, y = AG_sum.aq_agbds, by = c("marker", "chr", "pos", "ci.low", "ci.high", "lod", "pval", "pheno"), all =  T, sort = F)

chr.aq<-{}
pos.aq<-{}

#in below, replace "nrow" with "length"?
for(i in 1:nrow(sum.aq$lod)){
  chr.aq1<-sum.aq$lod[i,1]
  chr.aq2<-as.numeric(as.character(chr.aq1))
  chr.aq<-c(chr.aq, chr.aq1)
}
if(is.na(chr.aq)==TRUE) { chr.aq<-chr } else { chr.aq<-c(chr, chr.aq) }

for(i in 1:nrow(sum.aq$lod)){
  pos1.aq<-sum.aq$lod[i,2]
  pos.aq<-c(pos.aq, pos1.aq)
}
if(is.na(pos.aq)==TRUE) { pos.aq<-pos } else { pos.aq<-c(pos,pos.aq) }

qtl<-makeqtl(rils_84_impXO, chr=chr.aq, pos=pos.aq, what="prob")
createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
##########################################################################################

#(7) Uses forward selection and backward elimination model selection to probe the model space 
#for the best fit QTL model explaining your data.

##########################################################################################
#can calculate the penalities manually (but with scanTWO permutation results)
AG_scantwo.perms <- scantwo(AG_rils, pheno.col = 2:11, method = "em", model = "normal", n.perm =1000, verbose=TRUE, n.cluster = 10)
AG_calcpen = calc.penalties(AG_scantwo.perms)

AG_pen<-summary(AG_sub.perms) #pen should be a vector? 
#AG_pen are the 5% and 10% LOD thresholds; separate these by phenotype and use individually below?
AG.sw<-stepwiseqtl(AG_rils, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", penalties=AG_pen, model="normal", additive.only=TRUE)

rils_84_impXO.sw<-stepwiseqtl(rils_84_impXO, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", penalties=pen, model="normal", additive.only=TRUE)
swQTL<-capture.output(print(rils_84_impXO.sw))
swQTL2<-c(pheno, swQTL)
cat(swQTL2, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
cat(space, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
##########################################################################################

#(8) Adds additional QTL (if any) to the qtl object found by stepwise model selection.

##########################################################################################
sum.sw<-summary(rils_84_impXO.sw)
chr.sw<-{}
pos.sw<-{}
for(i in 1:nrow(sum.sw)){
  chr.sw1<-sum.sw[i,2]
  chr.sw2<-as.numeric(as.character(chr.sw1))
  chr.sw<-c(chr.sw, chr.sw2)
}
for(i in 1:nrow(sum.sw)){
  pos1.sw<-sum.sw[i,3]
  pos.sw<-c(pos.sw, pos1.sw)
}

qtl2<-makeqtl(rils_84_impXO, chr=chr.sw, pos=pos.sw, what="prob")
createqtl<- paste("Q", 1:qtl2$n.qtl, sep="")
### This formula with collapse "*" allows to check for interaction if there are two or more QTL for a trait
formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "*"))) ## changing "+" to "*"
rqtl<-refineqtl(rils_84_impXO, pheno.col=phenotype, qtl=qtl2, method="hk", model="normal")
##########################################################################################

#(9) Writes a file containing the peak marker and the flanking markers representing the 1.5 LOD interval of each QTL.
## i changed to 1.0 LOD drop
##########################################################################################
Q<-"Q"
space<-" "
for (i in 1:rqtl$n.qtl){
  interval<-capture.output(lodint(rqtl, qtl.index=i, drop=1.0, expandtomarkers=TRUE))
  q<-paste(Q, i, sep="")
  interval.new<-c(pheno, q, interval, space)
  cat(interval.new, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
}
##########################################################################################

#(10) Writes a csv file containing the results of ANOVA for the full and reduced models, 
#the % variance explained by each QTL and the estimated effect size 
#(half the distance between the means for each genotype class in the case of RILs).

#*NOTE*# You may receive a warning here about dropping individuals with missing phenotypes. 
#This is expected if such a case exists and does not effect the output of the code.

##########################################################################################
rils_84_impXO.ests<-fitqtl(rils_84_impXO, pheno.col=phenotype, qtl=rqtl, formula=formula, method="hk", dropone=TRUE, get.ests=TRUE, model="normal")
ests<-capture.output(summary(rils_84_impXO.ests))
ests<-c(pheno, ests)
write(ests, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
##########################################################################################

#(11) Generates a pdf file of the effect plots and a text file with means and standard error for each genotype class. 
#Generates a separate pdf file for every QTL.

#*NOTE*# This will likely produce a warning indicating that column names are being appended to file.

##########################################################################################
for(i in 1:length(chr.sw)) {
  b<-paste("Q",i, ".pdf",sep="")
  file<-paste("18_Marker effect plots", pheno, b)
  mar<-find.marker(cross2, chr=chr.sw[i], pos=pos.sw[i])
  pdf(file=file, width=11, height=8.5)
  plotPXG(rils_84_impXO, marker=mar, pheno.col=phenotype) ### set the color to black and check
  dev.off()
  phenoqtl<-paste(pheno, b)
  pheno.eff<-c(phenoqtl, space)
  means<-effectplot(cross2, pheno.col=phenotype, mname1=mar, draw=FALSE)
  cat(pheno.eff, file="19_means and SE.txt", sep="\n", append=TRUE)
  write.table(means$Means, file="19_means and SE.txt", sep=",", col.names="Means", row.names=TRUE, append=TRUE)
  cat(space, file="19_means and SE.txt", sep="\n", append=TRUE)
  write.table(means$SEs, file="19_means and SE.txt", sep=",", col.names="Standard Error", row.names=TRUE, append=TRUE)
  cat(space, file="19_means and SE.txt", sep="\n", append=TRUE)
}
##########################################################################################

#(12) Contingency statement if no QTL are detected in the initial genome scan for a given phenotype. You can ignore this step if you are not running the entire model selection section as a loop. 

##########################################################################################
{ else { null<-"	There were no LOD peaks above the threshold"
cat(pheno, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
cat(null, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
cat(space, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
cat(pheno, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
cat(null, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
cat(space, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
cat(pheno, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
cat(null, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
cat(space, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
cat(pheno, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
cat(null, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
cat(space, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
cat(pheno, file="19_means and SE.txt", sep="\n", append=TRUE)
cat(null, file="19_means and SE.txt", sep="\n", append=TRUE)
cat(space, file="19_means and SE.txt", sep="\n", append=TRUE) }
}
