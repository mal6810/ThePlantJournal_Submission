###QTL Figures #####
#################
### by JT Lovell
###A general pipeline for looking at mutli-trait multiple-qtl models
#http://rstudio-pubs-static.s3.amazonaws.com/223539_a394a8602be14af996e6f25b43c7f09e.html
#https://groups.google.com/forum/#!topic/Rqtl-disc/NjuBsMFXxeM

library(qtlTools)
#If getting an error when installing qtlTools due to update issue, 
#install it using the following commands:
library(devtools)
install_github("jtlovell/qtlTools")
library(qtlTools)

library(qtl)
library(ggplot2)
library(lsmeans)

phes_AG <- phenames(AG.rils) #collects phenotype names for ALL columns  
phes_AG <- phes_AG[c(2:11)] #collects phenotype names for selected columns (i.e. remove "id")

rils_84_impXO.sc1 <- scanone(rils_84_impXO, pheno.col = 1:33, method = "em", model = "normal")
#s1<-scanone(cross, pheno.col = phes, method = "hk")
cols<-rainbow(length(phes_AG)) ## this is for color pallete

cols<- c("red","blue","black") ## this is for color pallete

perms <- rils_84_impXO.perms

## to plot all QTLs for simple interval mapping (follow ALL steps)
plot(rils_84_impXO.sc1, type="n", ylim = c(0,max(as.matrix(rils_84_impXO.sc1[,-c(1:2)]))), 
     ylab = "LOD Score", main = "All phenotypes combined")
for(i in 1:length(phes)) plot(rils_84_impXO.sc1, add=T, lodcolumn = i, col = cols[i])
# to add threshold for all traits
add.threshold(rils_84_impXO.sc1, perms=perms[,c(1:33)], alpha=0.05, lodcolumn=i, gap=0, col="black", lty=2) # adds threshold line


## QTL peaks from scanones

perms <- rils_84_impXO.perms
perms<-scanone(rils_84_impXO, pheno.col = phes, method = "em", n.perm = 1000, verbose=F,n.cluster = 10)
sigqtl <- print(pullSigQTL(rils_84_impXO, pheno.col = phes,
                           s1.output = rils_84_impXO.sc1,  perm.output = perms, 
                           returnQTLModel=FALSE, alpha = 0.05,
                           controlAcrossCol = FALSE))

write.table(sigqtl, file = "scanone_simple_IM_significant_QTLs_alpha_0.05.txt", row.names = F, sep = "\t", quote = F)

## for confidence interval
cis<-calcCis(rils_84_impXO,s1.output=rils_84_impXO.sc1, perm.output=perms)
#this needs: cross file, scaneone file, and outout from scanone permutation
#lodint is default (the interval in which the LOD score is within some value of the QTL's maximum LOD score)
  #this differs for each QTL, as it is based on QTL-specifc parameters
head(cis)
write.table(cis, file = "scanone_simple_IM_QTLs.txt", row.names = F, sep = "\t", quote = F)
#provides high and low positions (in cM) for each QTL
#Leaf Temp QTL on Chr6 is just above the 0.05 alpha threshold (is 0.053)
#Adjust threshold in line 50 to include this QTL AND when calculating scanone.perm

## LOD threshold from permutation calculations..
sumx <- summary(perms[,c(1:2)]) #the 1:2 only gives LOD scores at 5% and 10% alpha for first two phenos
sumx = summary(perms) #gives LOD thresholds for all phenotypes in file
write.table(sumx, file = "scanone_simple_IM_QTLs_perm_intrvl.txt", row.names = T, sep = "\t", quote = F)


#### Write the genetic marker positions

gmpos <- as.data.frame(rils_84_impXO.sc1)

write.table(gmpos, file = "scanone_simple_IM_QTLs_genetic_map_positions.txt", row.names = T, sep = "\t", quote = F)


# to plot on genetic map
## margins, bottom,left, top, right
pdf(file="QTLPlot.pdf")
par(mar=c(5.1,4.1,4.1,4.1))
segmentsOnMap(rils_84_impXO, calcCisResults = cis, legendCex = .50, peaklod = NA, palette = rainbow, legendPosition = "bottomright",
              showPeaks = T, leg.inset = .00, tick.width = .05, chrBuffer = c(.05, .15),lwd = "byLod",
              leg.lwd = 2, max.lwd = 4, min.lwd = 1, title.adj =1.5, orderBy = "lod")
dev.off()

#adjustments for MAL figures
#change chrBuffer to 0.15, 0.35 (adds space between chromosomes and lines denoting QTL)
#use legendCex = 0.75 when saving as PDF
#palette = rainbow_hcl or qualitative_hcl (from the package "colorspace")
#legendPosition = "bottomright"
#lwd = "3"; avoid doing this "byLOD"
#max.lwd = 4
library(colorspace)
segmentsOnMap(AG_rils/RSA_rils, calcCisResults = cis, legendCex = .75, peaklod = NA, palette = rainbow_hcl, legendPosition = "bottomright",
              showPeaks = T, leg.inset = .00, tick.width = .05, chrBuffer = c(.15, .35),lwd = 3,
              leg.lwd = 4, max.lwd = 4, min.lwd = 1, title.adj =1.5, orderBy = "lod")

###using Haley Knot Method 
cross<-rils_84_impXO
cross<-calc.genoprob(cross)
phes <- phenames(rils_84_impXO) #same as above
phes <- phes[c(1:33)] #same as above
s1<-scanone(cross, pheno.col = phes, method = "hk")
cols<-rainbow(length(phes))
plot(s1, type="n", ylim = c(0,max(as.matrix(s1[,-c(1:2)]))), 
     ylab = "LOD Score", main = "all phenotypes plotted together")
for(i in 1:length(phes)) plot(s1, add=T, lodcolumn = i, col = cols[i])


perms<-scanone(cross, pheno.col = phes, method = "hk", n.perm = 1000, verbose=F,n.cluster = 10)
print(pullSigQTL(cross, pheno.col=phes,
                 s1.output = s1,  perm.output = perms, 
                 returnQTLModel=FALSE, alpha = 0.05,
                 controlAcrossCol = TRUE))
mods<-pullSigQTL(cross, pheno.col=phes,
                 s1.output = s1,  perm.output = perms, 
                 returnQTLModel=TRUE, alpha = 0.05,
                 controlAcrossCol = TRUE)
mods
cis<-calcCis(cross,s1.output=s1, perm.output=perms)
head(cis)
pdf(file="QTLPlot.pdf")
segmentsOnMap(cross, calcCisResults = cis, legendCex = .75,palette = rainbow_hcl, legendPosition = "bottomright",
              showPeaks = T, leg.inset = .00, tick.width = .05, chrBuffer = c(.15, .35),lwd = 2)
dev.off()
