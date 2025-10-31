library("ggplot2")
library("gridExtra")
library("tidyverse") # Package version: 1.2.1
library("ggthemes") # Package version: 4.0.1
library("scales") # Package version 1.0.0
library("RColorBrewer") # Package version 1.1-2
library("evaluate") # devtools::install_github('hadley/evaluate')
library("data.table") # Package version 1.12.0
library("rvest") # Package version 0.3.2
library("data.table")
library("limma")
library("yaml")

config <- read_yaml("penncnv_config.yaml")

LRR_SD.thres=config$thresholds[1]
NCNV.thres=config$thresholds[2]
WF.thres=config$thresholds[3] #Penn CNV suggest 0.04 might be a reasonable cutoff

# from Jacquemont paper log R ratio-standard deviation <0.35; B allele frequency-standard deviation <0.08; |waviness factor| <0.05, and CNVs <50 observations per array.

cnv.qual=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/FullDataTable_reclustered.qcsum",header=T,stringsAsFactors = F)
cnv.include=read.table("C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/FullDataTable_reclustered.qcpass",header=F)

#####################################################################################################
### We get a list of indiviuals who have failed and passed QC based on our predefined parameters
exclude.individuals=cnv.qual[which(as.numeric(as.character(cnv.qual$NumCNV)) >= NCNV.thres | as.numeric(as.character(cnv.qual$WF)) >=WF.thres | as.numeric(as.character(cnv.qual$WF))<=-WF.thres | as.numeric(as.character(cnv.qual$LRR_SD)) >=LRR_SD.thres),]
exclude.individuals$GROUP="QC Fail"
include.individuals=cnv.qual[-which(as.numeric(as.character(cnv.qual$NumCNV)) >= NCNV.thres | as.numeric(as.character(cnv.qual$WF)) >=WF.thres | as.numeric(as.character(cnv.qual$WF))<=-WF.thres | as.numeric(as.character(cnv.qual$LRR_SD)) >=LRR_SD.thres),]
include.individuals$GROUP="QC Pass"
cnv.qual=rbind(exclude.individuals,include.individuals)
#####################################################################################################
plot1=ggplot(cnv.qual,aes(x=WF,y=NumCNV,color=GROUP))+
  geom_point(shape=1)+
  xlab("Wavefactor")+
  ylab("N CNVs")+
  scale_color_manual(values=c("#E69F00","#56B4E9"))+
  theme(legend.position="none")

plot2=ggplot(cnv.qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
  geom_point(shape=1)+
  xlab("LRR SD")+
  ylab("N CNVs")+
  scale_color_manual(values=c("#E69F00","#56B4E9"))+
  theme(legend.position="bottom")


plot3=ggplot(cnv.qual, aes(x=BAF_drift)) + geom_histogram(alpha=.5, position="identity",colour="black", fill="white")+
  xlab("B-Allele Frequency Drift")+
  ylab("Frequency")

plot4=ggplot(cnv.qual, aes(x=LRR_SD)) + geom_histogram(colour="black", fill="white")+
  xlab("LRR SD")+
  ylab("Frequency")

no.outliers.removed=grid.arrange(arrangeGrob(plot1,plot2,ncol=1),arrangeGrob(plot3,plot4,ncol=1),ncol=2, widths=c(1,1))

#################################################################################

LRR_SD.vd=cnv.qual$LRR_SD>=as.numeric(LRR_SD.thres)
NCNV.vd=cnv.qual$NumCNV>=as.numeric(NCNV.thres)
WF.vd=(cnv.qual$WF<=(-1*as.numeric(WF.thres)) |  cnv.qual$WF>= as.numeric(WF.thres))
vd.all=cbind(LRR_SD.vd,NCNV.vd,WF.vd)
venn.diagram.lowest.thres=vennCounts(vd.all)
vennDiagram(vd.all,names = c(paste("LRR_SD"),paste("NCNV"),paste("WF")))
title(paste("Individuals excluded using the following parameters: \nLRR>",LRR_SD.thres," ; NCNVs>",NCNV.thres," ; WF>",WF.thres,sep=""))
#######################################################################################

exclude.individuals=cnv.qual[which(as.numeric(as.character(cnv.qual$NumCNV)) >= NCNV.thres | as.numeric(as.character(cnv.qual$WF)) >=WF.thres | as.numeric(as.character(cnv.qual$WF))<=-WF.thres | as.numeric(as.character(cnv.qual$LRR_SD)) >=LRR_SD.thres),]
cnv.qual=cnv.qual[! cnv.qual$File %in% unlist(exclude.individuals$File),]
if(length(unique(cnv.qual$GROUP))==1){
  plot1no=ggplot(cnv.qual,aes(x=WF,y=NumCNV,color=GROUP))+
    geom_point(shape=1)+
    xlab("Wavefactor")+
    ylab("N CNVs")+
    scale_color_manual(values=c("#56B4E9"))+
    theme(legend.position="none")
} else if(length(unique(cnv.qual$GROUP))==2){
  plot1no=ggplot(cnv.qual,aes(x=WF,y=NumCNV,color=GROUP))+
    geom_point(shape=1)+
    xlab("Wavefactor")+
    ylab("N CNVs")+
    scale_color_manual(values=c("#E69F00","#56B4E9"))+
    theme(legend.position="none")
}

if(length(unique(cnv.qual$GROUP))==1){
  plot2no=ggplot(cnv.qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
    geom_point(shape=1)+
    xlab("LRR SD")+
    ylab("N CNVs")+
    scale_color_manual(values=c("#56B4E9"))+
    theme(legend.position="bottom")
} else if(length(unique(cnv.qual$GROUP))==2){
  plot2no=ggplot(cnv.qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
    geom_point(shape=1)+
    xlab("LRR SD")+
    ylab("N CNVs")+
    scale_color_manual(values=c("#E69F00","#56B4E9"))+
    theme(legend.position="bottom")
}

plot3no=ggplot(cnv.qual, aes(x=BAF_drift)) + geom_histogram(alpha=.5, position="identity",colour="black", fill="white")+
  xlab("B-Allele Frequency Drift")+
  ylab("Frequency")
#geom_vline(xintercept=sd(cnv.qual[which(cnv.qual$GROUP=="PennCNV Pass"),]$BAF_drift)*2, linetype="dashed", color = "Red")

plot4no=ggplot(cnv.qual, aes(x=LRR_SD)) + geom_histogram(colour="black", fill="white")+
  xlab("LRR SD")+
  ylab("Frequency")

grid.arrange(arrangeGrob(plot1no,plot2no,ncol=1),arrangeGrob(plot3no,plot4no,ncol=1),ncol=2, widths=c(1,1))
#############################################

###added by Jo H to combine datafiles

folder_path <- "C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Genotyping/OutputfromPeter/FastFile-Q8AuViEDi3rmzjBc/FullDataTable_reclustered/FullDataTable_reclustered/output/qc_cnvs"

# Get a list of all file names in the folder (e.g., CSV files)
file_list <- list.files(path = folder_path, pattern = "\\.goodcnv$", full.names = TRUE)

# Use lapply to read all files with fread and bind them together
combined_goodcnv <- rbindlist(lapply(file_list, fread, header = FALSE, fill = TRUE))

###################################################################################################

### Sample /CNV exclusions

CNV.Calls=as.data.frame(combined_goodcnv)
names(CNV.Calls)=c("COORDS","NPROBES","SIZE","TYPE","ID","START_PROBE","END_PROBE","CONF")
CNV.Calls$NPROBES=gsub("numsnp=","",CNV.Calls$NPROBES)
CNV.Calls$SIZE=gsub("length=","",CNV.Calls$SIZE)
CNV.Calls$SIZE=gsub(",","",CNV.Calls$SIZE)
CNV.Calls$CHR=gsub("chr","",do.call(rbind,strsplit(CNV.Calls$COORDS,split=":"))[,1])
CNV.Calls$CONF=gsub("conf=","",CNV.Calls$CONF)
CNV.Calls$CONF=as.numeric(as.character(CNV.Calls$CONF))
CNV.Calls$START=as.numeric(gsub(",","",do.call(rbind,strsplit(gsub("chr","",do.call(rbind,strsplit(CNV.Calls$COORDS,split="-"))[,1]),split=":"))[,2]))
CNV.Calls$END=as.numeric(gsub(",","",do.call(rbind,strsplit(CNV.Calls$COORDS,split="-"))[,2]))
CNV.Calls=CNV.Calls[which(as.numeric(as.character(CNV.Calls$SIZE))>100000) & as.numeric(as.character(CNV.Calls$NPROBES))>=20 & CNV.Calls$CONF>=10,]
#CNV.Calls.Exclude=CNV.Calls[which(as.numeric(as.character(CNV.Calls$SIZE))<20000 | as.numeric(as.character(CNV.Calls$NPROBES))<10 | CNV.Calls$CONF<10),]
CNV.Calls.Exclude=CNV.Calls[which(as.numeric(as.character(CNV.Calls$SIZE))<100000 | as.numeric(as.character(CNV.Calls$NPROBES))<10 | CNV.Calls$CONF<10),]

#CNV.Calls=CNV.Calls[grep("_R1$",CNV.Calls$ID),]

CNV.Calls.HQ=CNV.Calls[!CNV.Calls$ID %in% unlist(exclude.individuals$File),] ### this is file with QC fails excluded I think
CNV.Calls.QCFail =CNV.Calls[CNV.Calls$ID %in% unlist(exclude.individuals$File),]##Jo H

hundredKB=CNV.Calls[which(as.numeric(CNV.Calls.HQ$SIZE) > 100000),]
hundredKB=CNV.Calls

########neuroCNVS##########################################################################



cnv.neurodev=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/CNVS.KK.2019.Sorted.txt",sep="\t",stringsAsFactors=F)
cnv.neurodev$CHR=gsub("chr","",do.call(rbind,strsplit(cnv.neurodev$V3,split=":"))[,1])
cnv.neurodev$START=as.numeric(gsub(",","",do.call(rbind,strsplit(gsub("chr","",do.call(rbind,strsplit(cnv.neurodev$V3,split="-"))[,1]),split=":"))[,2]))
cnv.neurodev$END=as.numeric(gsub(",","",do.call(rbind,strsplit(cnv.neurodev$V3,split="-"))[,2]))
cnv.neurodev$CLASS=do.call(rbind,lapply(1:nrow(cnv.neurodev),function(x){
  ifelse(length(grep("del",cnv.neurodev[x,1]))>0,yes=1,no=3)
}))
names(cnv.neurodev)[2]=c("LOCUS")
cnv.neurodev=cnv.neurodev[,c(1,4:7)]

##################################################################################

cnv.neuro.beta=unique(as.data.frame(do.call(rbind,lapply(1:nrow(CNV.Calls.HQ),function(x){
#cnv.neuro.beta=unique(as.data.frame(do.call(rbind,lapply(1:nrow(CNV.Calls.QCFail),function(x){
  cnv.chr=CNV.Calls.HQ[x,9]
  cnv.start=CNV.Calls.HQ[x,10]
  cnv.end=CNV.Calls.HQ[x,11]
  
  ### CNV spans entire region
  cnv.neurodev.match=cnv.neurodev[which(cnv.neurodev$CHR==cnv.chr & cnv.start<=cnv.neurodev$START & cnv.end>=cnv.neurodev$END),] 
  
  if(nrow(cnv.neurodev.match)>0){
    cnv.prop=1
    cnv.neurodev.match=cbind(cnv.neurodev.match,cnv.prop)
  } else if(nrow(cnv.neurodev.match)==0){
    
    ### CNV start is less than the start of the nd cnv and end of the CNV is less than than end of the nd cnv
    cnv.neurodev.match=cnv.neurodev[which(cnv.neurodev$CHR==cnv.chr & cnv.start<=cnv.neurodev$START & cnv.end<=cnv.neurodev$END & cnv.end>=cnv.neurodev$START),]
    
    if(nrow(cnv.neurodev.match)>0){
      cnv.neurodev.match.start=cnv.neurodev.match$START
      cnv.neurodev.match.end=cnv.neurodev.match$END
      
      cnv.prop=(cnv.end-cnv.neurodev.match$START)/(cnv.neurodev.match.end-cnv.neurodev.match.start)
      cnv.neurodev.match=cbind(cnv.neurodev.match,cnv.prop)
      ### CNV start is greater than start of nd cnv and ends after the end of the nd region
    } else if(nrow(cnv.neurodev.match)==0){
      cnv.neurodev.match=cnv.neurodev[which(cnv.neurodev$CHR==cnv.chr & cnv.start>=cnv.neurodev$START & cnv.end>=cnv.neurodev$END & cnv.start<=cnv.neurodev$END),]
      
      if(nrow(cnv.neurodev.match)>0){
        cnv.neurodev.match.start=cnv.neurodev.match$START
        cnv.neurodev.match.end=cnv.neurodev.match$END
        
        cnv.prop=(cnv.neurodev.match$END-cnv.start)/(cnv.neurodev.match.end-cnv.neurodev.match.start)
        cnv.neurodev.match=cbind(cnv.neurodev.match,cnv.prop)
        
      } else if(nrow(cnv.neurodev.match)==0){
        ### CNV falls inside nd region
        cnv.neurodev.match=cnv.neurodev[which(cnv.neurodev$CHR==cnv.chr & cnv.start>=cnv.neurodev$START & cnv.end<=cnv.neurodev$END),]
        
        if(nrow(cnv.neurodev.match)>0){
          cnv.neurodev.match.start=cnv.neurodev.match$START
          cnv.neurodev.match.end=cnv.neurodev.match$END
          
          cnv.prop=(cnv.end-cnv.start)/(cnv.neurodev.match.end-cnv.neurodev.match.start)
          cnv.neurodev.match=cbind(cnv.neurodev.match,cnv.prop)
        }}}}
  
  if(nrow(cnv.neurodev.match)>0){
    cbind(CNV.Calls.HQ[x,],cnv.neurodev.match)
  }
  
}))))

cnv.neuro.beta$TYPE=do.call(rbind,strsplit(cnv.neuro.beta$TYPE,split="="))[,2]
cnv.neuro.beta=cnv.neuro.beta[which(cnv.neuro.beta$TYPE==cnv.neuro.beta$CLASS),]

cnv.neuro.beta.fifpc=cnv.neuro.beta




###################################################################################################

#geneloc=read.table(file=paste("~/Hybrid.Magma/Data/NCBI37.3.gene.loc",sep=""),header=F,stringsAsFactors = F)
geneloc=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/NCBI37.3.gene.loc"
                     ,header=F,stringsAsFactors = F)


cnv.neuro.beta.genes=as.data.frame(do.call(rbind,lapply(1:nrow(cnv.neuro.beta.fifpc),function(x){
  a=cnv.neuro.beta.fifpc[x,]
  cnv.chr=a[1,9]
  cnv.start=a[1,10]
  cnv.end=a[1,11]
  
  geneloc.chr=geneloc[which(geneloc$V2==cnv.chr),]
  
  ### CNV covers whole gene
  gene.whole=geneloc.chr[which(cnv.start<geneloc.chr$V3 & cnv.end>geneloc.chr$V4),]
  gene.start=geneloc.chr[which(cnv.start<geneloc.chr$V3 & cnv.end<geneloc.chr$V4 & cnv.end>geneloc.chr$V3),]
  gene.end=geneloc.chr[which(cnv.start>geneloc.chr$V3 & cnv.end>geneloc.chr$V4 & cnv.start<geneloc.chr$V4),]
  gene.mid=geneloc.chr[which(cnv.start>geneloc.chr$V3 & cnv.end<geneloc.chr$V4),]
  
  
  if(nrow(gene.whole)>0 | nrow(gene.start)>0 | nrow(gene.end)>0 | nrow(gene.mid)>0){
    a$GENES=paste(unique(c(gene.whole$V6,gene.start$V6,gene.end$V6,gene.mid$V6)),collapse="|")
    a
  } else{
    a$GENES=""
    a
  }
})))

##########################################################

library(tidyverse)

Checks<- CNV.Calls.HQ%>%
  dplyr::filter(!ID %in% cnv.neuro.beta.genes$ID)

### Check for duplicate IDs. It would be unusual for an individual to have two different ND CNVs
### Most likely cause for dupliucates are having a smaller nested CNV within a larger one, where both are ND 
### Loop over the IDs, check to see if the CNV coordinates are duplicates, and if so just report the largest CNV. 

#if(nrow(cnv.neuro.beta.genes[which(duplicated(cnv.neuro.beta.genes$ID)),])>0){

#nested.cnv.dups=as.data.frame(do.call(rbind,lapply(cnv.neuro.beta.genes[which(duplicated(cnv.neuro.beta.genes$ID)),]$ID,function(x){
#  cnv.dups=cnv.neuro.beta.genes[which(cnv.neuro.beta.genes$ID==x),]

#  if(length(unique(cnv.dups$COORDS))==1){

### Get the smaller of the two CNVs (the smaller one will be the nested one)
#    cnv.nested=which.min(cnv.dups[,15]-cnv.dups[,14])
#    cnv.exclude=cnv.dups[cnv.nested,c(5,12)]
#  }
#})))
#}

#if(exists("nested.cnv.dups")==T){
#cnv.neuro.beta.fifpc=cnv.neuro.beta.genes[!(cnv.neuro.beta.genes$ID %in% unlist(nested.cnv.dups[,1]) & cnv.neuro.beta.genes$V1 %in% unlist(nested.cnv.dups[,2])),]
#}

#miscfolder="~/Downloads/"
#### NRXN1 exon processing
nrxn1=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/exon.NCBI.NRXN1",header=T,stringsAsFactors = F)
nrx1.exon.start=unlist(strsplit(nrxn1$exonStarts,split = ","))
nrx1.exon.end=unlist(strsplit(nrxn1$exonEnds,split = ","))
nrx1.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(nrx1.exon.start),function(x){
  nrx1.exon.start.loop=nrx1.exon.start[x]
  nrx1.exon.end.loop=nrx1.exon.end[x]
  cbind(nrx1.exon.start.loop,nrx1.exon.end.loop)
}))))
nrx1.exon.coords$EXON=paste("NRXN1 Exon ",seq(1:nrow(nrx1.exon.coords)),sep="")
###

### YWHAE exon processing #17p
YWHAE=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/exon.YWHAE",header=T,stringsAsFactors = F,fill = T)
YWHAE.exon.start=unlist(strsplit(YWHAE$exonStarts,split = ","))
YWHAE.exon.end=unlist(strsplit(YWHAE$exonEnds,split = ","))
YWHAE.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(YWHAE.exon.start),function(x){
  YWHAE.exon.start.loop=YWHAE.exon.start[x]
  YWHAE.exon.end.loop=YWHAE.exon.end[x]
  cbind(YWHAE.exon.start.loop,YWHAE.exon.end.loop)
}))))
YWHAE.exon.coords$EXON=paste("YWHAE Exon ",seq(1:nrow(YWHAE.exon.coords)),sep="")
###

### PAFAH1B1 exon processing #17p
PAFAH1B1=read.table(file="C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/exon.PAFAH1B1",header=T,stringsAsFactors = F,fill = T)
PAFAH1B1.exon.start=unlist(strsplit(PAFAH1B1$exonStarts,split = ","))
PAFAH1B1.exon.end=unlist(strsplit(PAFAH1B1$exonEnds,split = ","))
PAFAH1B1.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(PAFAH1B1.exon.start),function(x){
  PAFAH1B1.exon.start.loop=PAFAH1B1.exon.start[x]
  PAFAH1B1.exon.end.loop=PAFAH1B1.exon.end[x]
  cbind(PAFAH1B1.exon.start.loop,PAFAH1B1.exon.end.loop)
}))))
PAFAH1B1.exon.coords$EXON=paste("PAFAH1B1 Exon ",seq(1:nrow(PAFAH1B1.exon.coords)),sep="")
###

cnv.neuro.beta.genes$CRITERIA_MET=0
cnv.neuro.beta.genes$SIZE=as.numeric(cnv.neuro.beta.genes$SIZE)

cnvs.unique=unique(cnv.neuro.beta.genes$V1)

cnv.patho.criteria=as.data.frame(do.call(rbind,lapply(cnvs.unique,function(x){
  #cnv.patho.criteria[1,12]
  
  CNV.ND=cnv.neuro.beta.genes[which(cnv.neuro.beta.genes$V1==x),]
  
  as.data.frame(do.call(rbind,lapply(1:nrow(CNV.ND),function(y){
    CNV.ND.Line=CNV.ND[y,] 
    
    # 1q21.1	del/dup		 Size	>50%	of	critical	region
    if((x=="1q21.1 del" | x=="1q21.1 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 1p36 del/dup		 Size	>50%	of	critical	region,	affecting	GABRD
    if((x=="1p36 del (GABRD)" | x=="1p36 dup (GABRD)") & length(grep("GABRD",CNV.ND.Line$GENES))>0 & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # TAR del/dup		 Size	>50%	of	critical	region
    if((x=="TAR del" | x=="TAR dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # NRXN1 del, must include complete deletion of at least one exon
    if((x=="NRXN1 del") & length(grep("NRXN1",CNV.ND.Line$GENES))>0){
      
      
      nrxn1.exon.count=as.data.frame(do.call(rbind,lapply(1:nrow(nrx1.exon.coords),function(z){
        CNV.ND.Line.Start=CNV.ND.Line[1,10]
        CNV.ND.Line.End=CNV.ND.Line[1,11]
        
        exon.whole=nrx1.exon.coords[which(as.numeric(as.character(nrx1.exon.coords[z,1]))>CNV.ND.Line.Start & as.numeric(as.character(nrx1.exon.coords[z,2]))<CNV.ND.Line.End),]
        
        if(nrow(exon.whole)>0){
          exon.whole
        }
        
      })))
      
      if(nrow(nrxn1.exon.count)>0){
        CNV.ND.Line$CRITERIA_MET=1    
      }
    }
    
    # 2q11.2	del/dup	Size	>50%	of	critical	region,	affecting	both LMAN2L and	ARID5A
    if((x=="2q11.2 del" | x=="2q11.2 dup") & length(grep("LMAN2L",CNV.ND.Line$GENES))>0 & length(grep("ARID5A",CNV.ND.Line$GENES))>0 & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 2q13	del/dup		 Size	>50%	of	critical	region		
    if((x=="2q13 del" | x=="2q13 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 2q13	del/dup	(NPHP1)		 Size	>50%	of	critical	region,	affecting	NPHP1
    if((x=="2q13 del (NPHP1)" | x=="2q13 dup (NPHP1)") & length(grep("NPHP1",CNV.ND.Line$GENES))>0 & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 2q21.1	del/dup		 Size	>50%	of	critical	region		
    if((x=="2q21.1 del" | x=="2q21.1 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 2q37	del	(HDAC4)		 Size	>50%	of	critical	region,	affecting	HDAC4
    if((x=="2q37 del (HDAC4)" | x=="2q37 dup (HDAC4)") & CNV.ND.Line$cnv.prop>0.5 & length(grep("HDAC4",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 3q29	del/dup		 Size	>50%	of	critical	region		
    if((x=="3q29	del" | x=="3q29	dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # WolfâHirschhorn	del/dup		 Size	>50%	of	critical	region		
    if((x=="WolfâHirschhor del" | x=="WolfâHirschhor dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # Sotos	Syn/5q35	dup		 Size	>50%	of	critical	region		
    if(x=="Sotos Syn/5q35 dup" & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # Williams	Beuren	Syn	del/dup	Size	>50%	of	critical	region		
    if((x=="Williams	Beuren Syn del" | x=="Williams Beuren Syn dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    # 8p23.1	del/dup		 At	least	1Mbp	of	critical	region		
    if((x=="8p23.1 del" | x=="8p23.1 dup") & CNV.ND.Line$SIZE>1000000){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 9q34	del/dup	(EHMT1)		 At	least	1Mbp	CNVs,	including	EHMT1
    if(x=="9q34 del (EHMT1)" & CNV.ND.Line$SIZE>1000000 & length(grep("EHMT1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 10q11.21q11.23	del/dup		 Size	>50%	of	critical	region		
    if((x=="10q11.21q11.23 del" | x=="10q11.21q11.23 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 10q23	del/dup		 At	least	1Mbp,	including	NRG3 and	GRID1
    if((x=="10q23 del" | x=="10q23 dup") & CNV.ND.Line$SIZE>1000000 & length(grep("NRG3",CNV.ND.Line$GENES))>0 & length(grep("GRID1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # PotockiâShaffer	Syn	del/11p11.2	dup (EXT2), size	>50%	of	critical	region,	including	EXT2
    
    if(x=="Potocki-Shaffer syndrome del (EXT2)" & CNV.ND.Line$cnv.prop>0.5 & length(grep("EXT2",CNV.ND.Line$GENES))){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 15q11.2	del/dup		 Size	>50%	of	critical	region		
    if((x=="15q11.2 del" | x=="15q11.2 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 15q11.2 del/dup BP1-BP2		 Size	>50%	of	critical	region		
    if((x=="15q11.2 del BP1-BP2" | x=="15q11.2 dup BP1-BP2") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # PWS	del/dup		 Full	critical	region,	~4Mbp	
    if((x=="(PWS/AS del" | x=="(PWS/AS dup") & CNV.ND.Line$SIZE>4000000 & CNV.ND.Line$cnv.prop>0.8){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 15q13.3 del BP4-BP5		 Size	>50%	of	critical	region
    if(x=="15q13.3 del BP4-BP5" & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 15q24	del/dup		 At	least	1Mbp	between	the	AâE	intervals	
    if((x=="15q24 del" | x=="15q24 dup") & CNV.ND.Line$SIZE>1000000){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    
    
    # 15q25	del/dup		 At	least	1Mbp	between	the	AâD	intervals		
    if(x=="15q25 del" & CNV.ND.Line$SIZE>1000000){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 16p13.11	del/dup		 Size	>50%	of	critical	region
    if((x=="16p13.11 del" | x=="16p13.11 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 16p12.1	del		 Size	>50%	of	critical	region
    if(x=="16p12.1 del" & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    
    
    # 16p11.2	distal	del/distal	dup		 Size	>50%	of	critical	region		
    if((x=="16p11.2 distal del" | x=="16p11.2 distal dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 16p11.2	del/dup		 Size	>50%	of	critical	region
    if((x=="16p11.2 del" | x=="16p11.2 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    
    
    # 17p13.3	del/dup	(YWHAE)		 Exonic	deletions;	whole	gene	duplications
    
    if((x=="17p13.3 del (YWHAE)" | x=="17p13.3	dup (YWHAE") & CNV.ND.Line$cnv.prop>0.5 & length(grep("YWHAE",CNV.ND.Line$GENES))>0){
      
      YWHAE.exon.count=as.data.frame(do.call(rbind,lapply(1:nrow(YWHAE.exon.coords),function(z){
        CNV.ND.Line.Start=CNV.ND.Line[1,10]
        CNV.ND.Line.End=CNV.ND.Line[1,11]
        
        exon.whole=YWHAE.exon.coords[which(as.numeric(as.character(YWHAE.exon.coords[z,1]))>CNV.ND.Line.Start & as.numeric(as.character(YWHAE.exon.coords[z,2]))<CNV.ND.Line.End),]
        
        if(nrow(exon.whole)>0){
          exon.whole
        }
        
      })))
      
      if(nrow(YWHAE.exon.count)>0){
        CNV.ND.Line$CRITERIA_MET=1    
      }
    }
    
    
    # 17p13.3	del/dup	(PAFAH1B1) Exonic	deletions;	whole	gene	duplications
    if((x=="17p13.3 del (PAFAH1B1)" | x=="17p13.3	dup (PAFAH1B1") & CNV.ND.Line$cnv.prop>0.5 & length(grep("PAFAH1B1",CNV.ND.Line$GENES))>0){
      PAFAH1B1.exon.count=as.data.frame(do.call(rbind,lapply(1:nrow(PAFAH1B1.exon.coords),function(z){
        CNV.ND.Line.Start=CNV.ND.Line[1,10]
        CNV.ND.Line.End=CNV.ND.Line[1,11]
        
        exon.whole=PAFAH1B1.exon.coords[which(as.numeric(as.character(PAFAH1B1.exon.coords[z,1]))>CNV.ND.Line.Start & as.numeric(as.character(PAFAH1B1.exon.coords[z,2]))<CNV.ND.Line.End),]
        
        if(nrow(exon.whole)>0){
          exon.whole
        }
        
      })))
      
      if(nrow(PAFAH1B1.exon.count)>0){
        CNV.ND.Line$CRITERIA_MET=1    
      }
    }
    # 17q11.2	del/dup	(NF1)	>50%	of	critical	region,	affecting	NF1
    if((x=="17q11.2 del (NF1)" | x=="17q11.2 dup (NF1)") & CNV.ND.Line$cnv.prop>0.5 & length(grep("NF1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    
    # SmithâMagenis/PotockiâLupski	Syndrome	Size	>50%	of	critical	region		
    
    if((x=="Smith-Magenis syndrome del" | x=="Potocki-Lupski syndrome dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 17q11.2	del/dup	(NF1)		 Size	>50%	of	critical	region,	affecting	NF1
    if((x=="17q11.2 del (NF1)" | x=="17q11.2 dup (NF1)") & CNV.ND.Line$cnv.prop>0.5 & length(grep("NF1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 17q12	del/dup		 Size	>50%	of	critical	region		
    if((x=="17q12 del" | x=="17q12 dup" | x=="Renal cysts and diabetes syndrome del") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 17q21.31	del/dup		 Size	>50%	of	critical	region		
    if((x=="17q21.31 del" | x=="17q21.31 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 17q23.1q23.2	del		 Size	>50%	of	critical	region		
    if(x=="17q23.1q23.2 del" & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 17p12       del              Size   >50%    of      critical        region
    if(x=="17p12 del" & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1
    }
    
    # 22q11.2	del/dup		 Size	>50%	of	critical	region		
    if((x=="22q11.2 del" | x=="22q11.2 dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # 22q11.2	distal	del/dup		 Size	>50%	of	critical	region
    if((x=="22q11.2 distal del" | x=="22q11.2 distal dup") & CNV.ND.Line$cnv.prop>0.5){
      CNV.ND.Line$CRITERIA_MET=1    
    }
    
    # SHANK3 del/dup	 At	least	1Mbp	CNVs,	including	SHANK3
    if((x=="SHANK3 del" | x=="SHANK3 dup") & length(grep("SHANK3",CNV.ND.Line$GENES))>0 & CNV.ND.Line$SIZE>1000000){
      CNV.ND.Line$CRITERIA_MET=1    
    } 
    
    
    CNV.ND.Line
    
  })))
})))

patho.criteria.met=cnv.patho.criteria[which(cnv.patho.criteria$CRITERIA_MET==1),]

### Deal with smaller nested CNVs 

patho.criteria.met.no.nested=as.data.frame(do.call(rbind,lapply(unique(patho.criteria.met$ID),function(x){
  
  ### Select row with individual and criteria met
  
  ### Selects individuals who have criteria =1 
  
  a=patho.criteria.met[which(patho.criteria.met$ID==x),]
  
  
  
  if(nrow(a)==1){
    
    a
    
  } else if(nrow(a)>1){
    
    
    
    ### Check nested TAR / 1q21, remove TAR from the results row
    
    if((length(grep("TAR",a$V1))+length(grep("1q21",a$V1)))==2){
      
      b1=a[-grep("TAR",a$V1),]
      
      b1
      
    }
    
    
    
    ### Check nested 15q11.2 in PWS/AS 
    
    
    
    if((length(grep("15q",a$V1))+length(grep("PWS",a$V1)))==2){
      
      b2=a[-grep("15q",a$V1),]
      
      b2
      
    }
    
    
    
    ### Check nested 16p11.2 distal within 
    
    
    
    if((length(grep("distal",a$V1))+length(grep("16p11.2 del",a$V1)))==2){
      
      b3=a[-grep("distal",a$V1),]
      
      b3
      
    }
    
    
    
    ###  A final statement to deal with what happens if there are still two separate records for an individual
    
    ### This just prints all lines. Note this does not deal with multiple ND CNVs at the same locus.
    
    if(exists("b1")==F & exists("b2")==F & exists("b3")==F ){
      
      a
      
    }  else if(exists("b1")==T){
      
      b1
      
    }
    
    else if(exists("b2")==T){
      
      b2
      
    } 
    
    else if(exists("b3")==T){
      
      b3
      
    } 
    
  }
  
})))

###

patho.criteria.met=patho.criteria.met.no.nested
patho.criteria.met=patho.criteria.met[order(as.numeric(as.character(patho.criteria.met[,9])),as.numeric(as.character(patho.criteria.met[,10]))),]
cnv.patho=as.data.frame(table(patho.criteria.met$V1))
cnv.patho=unique(merge(patho.criteria.met[,c(12,9,10)],cnv.patho,by.x="V1",by.y="Var1",sort=F)[,c(1,4)])


cnv.patho$DATASET="July25"
names(cnv.patho)=c("Neurodevelopmental CNV","N","DATASET")

write.table(rbind(patho.criteria.met.no.nested,cnv.patho.criteria[which(cnv.patho.criteria$CRITERIA_MET==0),]),file="C:\\Users\\sapjeh\\OneDrive - Cardiff University\\Sleep Detectives\\Family Environment Analysis\\Genotyping\\Pipeline\\July25.Patho.Annotated_2025-10-20.txt",col.names=T,row.names=F,quote=F,sep="\t")

CalledCNVS<- rbind(patho.criteria.met.no.nested,cnv.patho.criteria[which(cnv.patho.criteria$CRITERIA_MET==0),])

write.csv(CalledCNVS, ,file="C:\\Users\\sapjeh\\OneDrive - Cardiff University\\Sleep Detectives\\Family Environment Analysis\\Genotyping\\Pipeline\\July25.CNVCalls_2025-10-20.csv")
################################################################################################



invisible(lapply(unique(patho.criteria.met$ID),function(x){
  #browser()
  patho.id=patho.criteria.met[which(patho.criteria.met$ID==x),]
  
  if(nrow(patho.id)==1){
    patho.cnv.lower.coords=as.numeric(as.character(patho.id[,10]))-(patho.id[1,3]*1.3)
    patho.cnv.upper.coords=as.numeric(as.character(patho.id[,11]))+(patho.id[1,3]*1.3)
    patho.cnv.chr=as.numeric(as.character(patho.id[,9]))
    cnv.raw=as.data.frame(fread(input=paste("C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/non-iPS/",x,sep=""),header=T,sep="\t"))# all need to be in same folder
    cnv.raw.params=cnv.raw[which(cnv.raw$Position>=patho.cnv.lower.coords & cnv.raw$Position<=patho.cnv.upper.coords & cnv.raw$Chr==patho.cnv.chr),]
    cnv.raw.params$GROUP=c("Probes outside called CNV")
    cnv.raw.params[which(cnv.raw.params$Position>=as.numeric(as.character(patho.id[,10])) & cnv.raw.params$Position<=as.numeric(as.character(patho.id[,11])) & cnv.raw.params$Chr==patho.cnv.chr),]$GROUP=c("Probes within individually called CNV")
    cnv.raw.params$GROUP <- factor(cnv.raw.params$GROUP, levels = c("Probes within individually called CNV","Probes outside called CNV"))
    
    
    baf.col=grep("Allele",names(cnv.raw.params))
    logr.col=grep("Log",names(cnv.raw.params))
    
    baf=ggplot(data=cnv.raw.params,aes(x=cnv.raw.params[,3],y=cnv.raw.params[,baf.col],color=as.factor(GROUP))) +
      geom_point(shape=1) +
      xlab("Base Position")+
      ylab("B-Allele Frequency")+
      ggtitle(paste(x,patho.id$V1,sep=" "))+
      geom_hline(yintercept=0.5, linetype="dashed", color = "Black")+
      geom_vline(xintercept=patho.id[,14], linetype="dashed", color = "Green")+
      geom_vline(xintercept=patho.id[,15], linetype="dashed", color = "Green")+
      ylim(c(0,1))+
      labs(color="CNV Probe Legend")
    
    
    
    lrr=ggplot(data=cnv.raw.params,aes(x=cnv.raw.params[,3],y=cnv.raw.params[,logr.col],color=as.factor(GROUP))) +
      geom_point(shape=1) +
      xlab("Base Position")+
      ylab("Log R Ratio")+
      ylim(c(-1,1))+
      geom_hline(yintercept=0, linetype="dashed", color = "Black")+
      geom_vline(xintercept=patho.id[,14], linetype="dashed", color = "Green")+
      geom_vline(xintercept=patho.id[,15], linetype="dashed", color = "Green")+
      labs(color="CNV Probe Legend")
    
    
    
    grid.arrange(arrangeGrob(baf,lrr,ncol=1))
    
    png(paste("C:\\Users\\sapjeh\\OneDrive - Cardiff University\\Sleep Detectives\\Family Environment Analysis\\Genotyping\\Pipeline\\NeuroDevelopmentalPlots\\",x,".",patho.id$V1,".png",sep=""), width = 10, height = 4, units = 'in', res = 300)
    grid.arrange(arrangeGrob(baf,lrr,ncol=1))
    dev.off()
  } 
  
  else if(nrow(patho.id)>1){
    lapply(1:nrow(patho.id),function(y){
      
      patho.cnv.lower.coords=as.numeric(as.character(patho.id[y,10]))-(patho.id[y,3]*1.3)
      patho.cnv.upper.coords=as.numeric(as.character(patho.id[y,11]))+(patho.id[y,3]*1.3)
      patho.cnv.chr=as.numeric(as.character(patho.id[y,9]))
      cnv.raw=as.data.frame(fread(paste0("C:/Users/sapjeh/OneDrive - Cardiff University/Sleep Detectives/Family Environment Analysis/Genotyping/Pipeline/non-iPS/",x,sep=""),header=T,sep="\t"))
      cnv.raw.params=cnv.raw[which(cnv.raw$Position>=patho.cnv.lower.coords & cnv.raw$Position<=patho.cnv.upper.coords & cnv.raw$Chr==patho.cnv.chr),]
      cnv.raw.params$GROUP=c("Probes outside called CNV")
      cnv.raw.params[which(cnv.raw.params$Position>=as.numeric(as.character(patho.id[y,10])) & cnv.raw.params$Position<=as.numeric(as.character(patho.id[y,11])) & cnv.raw.params$Chr==patho.cnv.chr),]$GROUP=c("Probes within individually called CNV")
      cnv.raw.params$GROUP <- factor(cnv.raw.params$GROUP, levels = c("Probes within individually called CNV","Probes outside called CNV"))
      
      
      baf.col=grep("Allele",names(cnv.raw.params))
      logr.col=grep("Log",names(cnv.raw.params))
      
      baf=ggplot(data=cnv.raw.params,aes(x=cnv.raw.params[,3],y=cnv.raw.params[,baf.col],color=as.factor(GROUP))) +
        geom_point(shape=1) +
        xlab("Base Position")+
        ylab("B-Allele Frequency")+
        ggtitle(paste(x,patho.id$V1,sep=" "))+
        geom_hline(yintercept=0.5, linetype="dashed", color = "Black")+
        geom_vline(xintercept=patho.id[y,14], linetype="dashed", color = "Green")+
        geom_vline(xintercept=patho.id[y,15], linetype="dashed", color = "Green")+
        ylim(c(0,1))+
        labs(color="CNV Probe Legend")
      
      lrr=ggplot(data=cnv.raw.params,aes(x=cnv.raw.params[,3],y=cnv.raw.params[,logr.col],color=as.factor(GROUP))) +
        geom_point(shape=1) +
        xlab("Base Position")+
        ylab("Log R Ratio")+
        ylim(c(-1,1))+
        geom_hline(yintercept=0, linetype="dashed", color = "Black")+
        geom_vline(xintercept=patho.id[y,14], linetype="dashed", color = "Green")+
        geom_vline(xintercept=patho.id[y,15], linetype="dashed", color = "Green")+
        labs(color="CNV Probe Legend")
      
      grid.arrange(arrangeGrob(baf,lrr,ncol=1))
      
      png(paste("C:\\Users\\sapjeh\\OneDrive - Cardiff University\\Sleep Detectives\\Family Environment Analysis\\Genotyping\\Pipeline\\NeuroDevelopmentalPlots\\",x,".",patho.id$V1,".png",sep=""), width = 10, height = 4, units = 'in', res = 300)
      grid.arrange(arrangeGrob(baf,lrr,ncol=1))
      dev.off()
    })  
  }
  
}))
