##################################################################
##################################################################
## Processing and exploration of DArTseq data
##################################################################
##################################################################
## Requires *SNP_2.csv file supplied by DArT and matching sample
## *metadata.csv in a folder specific to the DArT run (within 
## the data directory). 
##
##################################################################
## Load packages
##################################################################

library(dartR)
library(StAMPP)
library(gplots) # only for heatmap.2?
# library(ggplot2)
library(plotly)
library(hierfstat)
library(poppr)

##################################################################
## Setup
##################################################################
## Specify main parameters
setwd(".")
overwrite=TRUE
## Directories: choose an indir and an analysis name
indir="data/OrderAppendix_1_DPro19-4334" ## Psubset
# indir="data/OrderAppendix_2_DPro19-4334" ## PgOnly
indir="data/OrderAppendix_3_DPro19-4334" ## PgRels
# analysis_name="Psubset"
# analysis_name="Pg"
analysis_name="PgRels"

## Calculate other parameters
dartname=list.files(indir,"SNP_2.csv")
dartfile=file.path(indir,dartname)
metadataname=list.files(indir,"metadata.csv")
metadatafile=file.path(indir,metadataname)
outdir=file.path("output",analysis_name)
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("temp")) dir.create("temp")
if(dir.exists(outdir)){
  if(overwrite){
    print("Overwriting existing directory")
  }else print("Directory already exists")
} else dir.create(outdir)

##################################################################
## Load data
##################################################################
## Read data file. Note that metadata are not read in here, as
##  the order must be the same, but the dartfile is not in a
##  sensible order.
mygl=gl.read.dart.2row(datafile = dartfile,
                       topskip = 6,nmetavar = 26)

## Read metadata file and sort according to sample order in dartfile
metadata0=read.csv(metadatafile, stringsAsFactors = F)
metadata=metadata0[match(mygl@ind.names,metadata0$id),]

## Add sample metadata to mygl
mygl@pop=factor(metadata$pop)
# mygl@pop=factor(metadata$pop,levels=c("gilesii","phylicifolia","monticola","Evans Crown"))


##################################################################
## Filter data
##################################################################
## Explore data quality
gl.report.bases(mygl)
gl.report.callrate(mygl)
gl.report.callrate(mygl,method = "ind")
# gl.report.hamming(mygl) # requires too much memory
gl.report.repavg(mygl)
gl.report.secondaries(mygl)
gl.report.monomorphs(mygl)
# gl.report.hwe(mygl) # produces a large dataframe

## Filter data - using the same object to limit memory usage
mygl1=gl.filter.repavg(mygl,threshold = 0.95)
nLoc(mygl1);nInd(mygl1)
mygl1=gl.filter.callrate(mygl1,method = "loc",threshold = 0.95)
mygl1=gl.filter.callrate(mygl1,method = "ind",threshold = 0.9)
nLoc(mygl1);nInd(mygl1)
mygl1=gl.filter.secondaries(mygl1)
nLoc(mygl1);nInd(mygl1)

mygl1=gl.recalc.metrics(mygl1)
pop0=mygl1@pop
save(mygl1,file = file.path("temp",paste0(analysis_name,"_filtered.rda")))

##################################################################
## Population-level statistics
##################################################################

## Compare heterozygosity among populations
ho_pop <- gl.report.heterozygosity(mygl1)
pdf(file.path(outdir,paste0(analysis_name,"_ho_pop.pdf")));gl.report.heterozygosity(mygl1);dev.off()
write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_pop.txt")))
# gl.report.maf(mygl) # takes a lot of time

## Examine private alleles by pairs of populations
pa_pw_pop <- gl.report.pa.pop(mygl1)
pa_pw_pop
write.table(pa_pw_pop,file = file.path(outdir,paste0(analysis_name,"_pa_pw_pop.txt")))

#hwe <- gl.report.hwe(mygl1,subset=c("gilesii","phylicifolia")) # 
# not sure how this function filters rows

## Genetic distances
mydist1 <- gl.dist.pop(mygl1,method = "pcfixed"); mydist1
mydist2 <- gl.dist.pop(mygl1,method = "pa"); mydist2
mydist3 <- gl.dist.pop(mygl1,method = "gower",upper = T); mydist3
mydist4 <- gl.dist.pop(mygl1,method = "euclidean",upper = T); mydist4
write.table(mydist1,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pcfixed.txt")),sep="\t")
write.table(mydist2,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pa.txt")),sep="\t")
write.table(as.matrix(mydist3),file = file.path(outdir,paste0(analysis_name,"_pop_dist_gower.txt")),sep="\t")
write.table(as.matrix(mydist4),file = file.path(outdir,paste0(analysis_name,"_pop_dist_euclidean.txt")),sep="\t")

## Nei's (1972) genetic distance

mygeno <- stamppConvert(mygl1,type = "genlight")
mydist5 <- stamppNeisD(mygeno,pop=T)
write.table(mydist5,file = file.path(outdir,paste0(analysis_name,"_pop_dist_Nei1972.txt")),sep="\t")

## F statistics with STAMPP
pwfst <-stamppFst(mygl1, nboots=1000, percent=95, nclusters=1)
pwfst$Fsts
pwfst$Pvalues
write.table(pwfst$Fsts,file = file.path(outdir,paste0(analysis_name,"_pop_fst.txt")),sep="\t")

## Heterozygosity and F statistics with hierfstat
sumstats=basic.stats(cbind(mygl1@pop,as.data.frame(mygl1)))
sumstats$overall
pop_summary=data.frame(N=colMeans(sumstats$n.ind.samp,na.rm=T),
                       Ho=colMeans(sumstats$Ho,na.rm=T),
                       He=colMeans(sumstats$Hs,na.rm=T),
                       Fis=colMeans(sumstats$Fis,na.rm=T))
write.table(pop_summary,file = file.path(outdir,paste0(analysis_name,"_pop_summary.txt")),sep="\t")

hist(sumstats$perloc$Ho)
plot(sumstats$perloc$Hs,sumstats$perloc$Ho)

if(analysis_name=="PgRels"){ ## Not working yet
  phyl_sumstats=basic.stats(cbind(mygl1@pop,as.data.frame(mygl1))[pop0=="phylicifolia",])
  gilesii_sumstats=basic.stats(cbind(mygl1@pop,as.data.frame(mygl1))[pop0=="gilesii",])
  evans_sumstats=basic.stats(cbind(mygl1@pop,as.data.frame(mygl1))[pop0=="Evans Crown",])
}

##################################################################
## Ordination
##################################################################
## Principal coordinates analysis
mypc <- gl.pcoa(mygl1,nfactors = 6)
## Scree plot
barplot(mypc$eig/sum(mypc$eig)*100)

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=2)
plot1=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=2)
plot2=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=3)
plot3=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=4)
plot4=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=5)
plot5=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=6)
plot6=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=2,yaxis=3)
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_12.pdf")));plot1;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_13.pdf")));plot2;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_14.pdf")));plot3;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_15.pdf")));plot4;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_16.pdf")));plot5;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_23.pdf")));plot6;dev.off()

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=2)
ggplotly(tooltip=c("ind","x","y")) # shows pop twice or not at all

## Doesn't work (yet):
plot1=gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=2)+
  geom_text(data = subset(data.frame(id=indNames(mygl1),mypc$scores),id=="TCW 577"),aes(label=id))


## TCW 557 is the odd one out in Prostanthera gilesii


##################################################################
## Individual-level statistics - in progress
##################################################################

## Genomic Relationship Calculation - not quite sure what to make of these
gmatrix=stamppGmatrix(mygeno)
gmatrix[1:10,1:10]
## Heat maps - not sure if these make sense
heatmap(gmatrix,symm=T)
heatmap.2(gmatrix,density="none",trace="none")

## Correlation matrix
cormat=cor(as.data.frame(mygl1),use="pairwise.complete.obs")


