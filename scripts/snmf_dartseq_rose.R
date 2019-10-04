##################################################################
##################################################################
## Perform SNMF on filtered DArTseq data
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
library(adegenet)
library(LEA)
library(poppr)


##################################################################
## Setup
##################################################################
## Specify main parameters
setwd(".")
overwrite=FALSE
## Directories: choose an indir and an analysis name
# indir="temp" 
# analysis_name="Psubset"
# analysis_name="Pg"
analysis_name="PgRels"

## Calculate other parameters
# dartname=list.files(indir,"SNP_2.csv")
# dartfile=file.path(indir,dartname)
# metadataname=list.files(indir,"metadata.csv")
# metadatafile=file.path(indir,metadataname)
# outdir=file.path("output",analysis_name)
# if(!dir.exists("output")) dir.create("output")
# if(!dir.exists("temp")) dir.create("temp")
# if(dir.exists(outdir)){
#   if(overwrite){
#     print("Overwriting existing directory")
#   }else print("Directory already exists")
# } else dir.create(outdir)

##################################################################
## Load data and convert to geno object
##################################################################
## Load the filtered data file, a genlight object called mygl1,
##  saved as a .rda file.
load(file = file.path("temp",paste0(analysis_name,"_filtered.rda")))
indNames(mygl1) <- gsub(" ","",indNames(mygl1))

## Filter and thin
# mygl2=popsub(mygl1,blacklist = "Walls")
# nLoc(mygl2);nInd(mygl2)
mygl2 <- mygl1
mygl2=gl.filter.repavg(mygl2,threshold = 0.95)
nLoc(mygl2);nInd(mygl2)
mygl2=gl.filter.callrate(mygl2,method = "loc",threshold = 0.95)
mygl2=gl.filter.callrate(mygl2,method = "ind",threshold = 0.9)
nLoc(mygl2);nInd(mygl2)
mygl2=gl.filter.secondaries(mygl2)
nLoc(mygl2);nInd(mygl2)
mygl2=gl.filter.monomorphs(mygl2)
nLoc(mygl2);nInd(mygl2)
mygl2=gl.filter.maf(mygl2, threshold = 0.05)
nLoc(mygl2);nInd(mygl2)

# mygl2=gl.recalc.metrics(mygl2)
mygl2=mygl2[order(mygl2@pop),]
pop0=mygl2@pop
save(mygl2,file = file.path("temp",paste0(analysis_name,"_filtered_maf.rda")))

## Write to STRUCTURE format
gl2structure(mygl2, 
             indNames = mygl2@ind.names,
             addcolumns = as.numeric(pop0),
             ploidy = 2,
             exportMarkerNames = FALSE,
             outfile = file.path("output",analysis_name,paste0(analysis_name,".struc")),
             # outpath = file.path("output",analysis_name),
             outpath = ".",
             v = 1)


## Read STRUCTURE file into genotype file 
mygeno=struct2geno(file.path("output",analysis_name,paste0(analysis_name,".struc")),
                   ploidy = 2,FORMAT = 2,
                   # extra.row = 1,
                   extra.column = 2)



##################################################################
## Find natural genetic distance thresholds
##################################################################
## First plot genetic distance histogram, and effects of different
##  thresholds and algorithms - needed to tweak input to get decent plot
hist(mlg.filter(mysc,threshold = 1e-6,stats = "DISTANCES"))
plot_filter_stats(mysc,filter_stats(mysc,threshold = .1),
                  mlg.filter(mysc,threshold=0,stats = "DISTANCES"),
                  breaks = "scott")
hist(mlg.filter(mysc,threshold = 1e-6,stats = "DISTANCES"),breaks = 100)
range(mlg.filter(mysc,threshold = 1e-6,stats = "DISTANCES"))

## Then see if the automated approach gives a similar result
cutoff_predictor(filter_stats(mysc)$farthest$THRESHOLDS)

## Cutoffs to use: 
auto_threshold=0.013
empir1_threshold=0.07
empir2_threshold=0.05
empir3_threshold=0.03
empir4_threshold=0.022

##################################################################
## Apply thresholds to data 
##################################################################
auto_mlg <- mlg.filter(mysc,threshold = auto_threshold)
auto_mysc <- mysc
auto_mysc$mlg <- auto_mlg
empir1_mlg <- mlg.filter(mysc,threshold = empir1_threshold)
empir1_mysc <- mysc
empir1_mysc$mlg <- empir1_mlg
empir2_mlg <- mlg.filter(mysc,threshold = empir2_threshold)
empir2_mysc <- mysc
empir2_mysc$mlg <- empir2_mlg
empir3_mlg <- mlg.filter(mysc,threshold = empir3_threshold)
empir3_mysc <- mysc
empir3_mysc$mlg <- empir3_mlg
empir4_mlg <- mlg.filter(mysc,threshold = empir4_threshold)
empir4_mysc <- mysc
empir4_mysc$mlg <- empir4_mlg



sclist <- list(auto_mysc,empir1_mysc,empir2_mysc,empir3_mysc,empir4_mysc)
mlglist <- sapply(sclist,mlg)
mlglist

##################################################################
## Explore clonal lineages
##################################################################

## plot on PCoA
mygl2@pop=as.factor(empir2_mlg)
mypc <- gl.pcoa(mygl2,nfactors = 6)
## Scree plot
barplot(mypc$eig/sum(mypc$eig)*100)
plot1=gl.pcoa.plot(mypc,mygl2,labels="pop",xaxis=1,yaxis=2)
plot1


## Plot phylo
plot.phylo(upgma(mlg.filter(mysc,threshold=0,stats = "DISTANCES")))

## Minimum spanning networks, with various thresholds
poppr.msn(auto_mysc,mlg.filter(mysc,threshold=0,stats = "DISTANCES"))
poppr.msn(empir1_mysc,mlg.filter(mysc,threshold=0,stats = "DISTANCES"))
poppr.msn(empir2_mysc,mlg.filter(mysc,threshold=0,stats = "DISTANCES"))
poppr.msn(empir3_mysc,mlg.filter(mysc,threshold=0,stats = "DISTANCES"))
poppr.msn(empir4_mysc,mlg.filter(mysc,threshold=0,stats = "DISTANCES"))

##################################################################
## Clonal diversity statistics
##################################################################

## Tabulate data
mytable=mlg.table(empir2_mysc)
mytable
diversity_stats(mytable)

## Function for clonal fraction (from https://cran.r-project.org/web/packages/poppr/vignettes/mlg.html)
myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x > 0)/rowSums(x)
  } else {                 # if it's a vector
    res <- sum(x > 0)/sum(x)
  }
  return(res)
}

## Table with clonal fraction added
diversity_stats(mytable, CF = myCF)

mytableCI <- diversity_ci(mytable, n = 100L, rarefy = FALSE, raw = FALSE, CF = myCF)

## Summary statistics
sumstats=basic.stats(cbind(empir2_mysc@pop,as.data.frame(empir2_mysc)))
sumstats$overall
pop_summary=data.frame(N=colMeans(sumstats$n.ind.samp,na.rm=T),
                       Ho=colMeans(sumstats$Ho,na.rm=T),
                       He=colMeans(sumstats$Hs,na.rm=T),
                       Fis=colMeans(sumstats$Fis,na.rm=T))
pop_summary

hist(sumstats$perloc$Ho)
plot(sumstats$perloc$Hs,sumstats$perloc$Ho)

sumstats0=basic.stats(cbind(mysc@pop,as.data.frame(mysc)))
sumstats0$overall
pop_summary0=data.frame(N=colMeans(sumstats0$n.ind.samp,na.rm=T),
                        Ho=colMeans(sumstats0$Ho,na.rm=T),
                        He=colMeans(sumstats0$Hs,na.rm=T),
                        Fis=colMeans(sumstats0$Fis,na.rm=T))
pop_summary0

plot(sumstats0$perloc$Hs,sumstats0$perloc$Ho)

##################################################################
## 
##################################################################

# poppr.amova(mysc,hier = ~pop)





