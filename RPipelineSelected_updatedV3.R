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

# The vast majority of this script is written by
#- "Leon Hubbard"
#- "Kimberley Kendall"
#- "Elliott Rees"
#- "George Kirov"
# Changes to make directories more accessable and modifications by Jo Haddon
# Yaml adaptations, figure updates Peter Richardson

#setwd("./CNV_data/CNV_repo/") # Just temporary while I work on getting this running. Then I will just call it from the current dir anyway.

# All of the variables to set here including addresses are set in the yaml file.
config <- read_yaml("penncnv_config.yaml")

LRR_SD_thres=config$thresholds[[1]]
NCNV_thres=config$thresholds[[2]]
WF_thres=config$thresholds[[3]] #Penn CNV suggest 0.04 might be a reasonable cutoff


# from Jacquemont paper log R ratio-standard deviation <0.35; B allele frequency-standard deviation <0.08; |waviness factor| <0.05, and CNVs <50 observations per array.
# Need to name for either HPC or local?
prefix <- config$thresholds$sample_prefix
cnv_qc_sum <- paste0(prefix, ".qcsum")
cnv_qc_include <- paste0(prefix, ".qcpass")

# TODO Need to update this to make it adjust accordingly
FOLDER_PATH <- paste0("./tmp_scratch/", prefix, "/output/qc_cnvs/")


cnv_kk_file <- "./penn_cnv_files/CNVS.KK.2019.Sorted.txt"


cnv_qual <- read.table(
  file = file.path(
    "./tmp_scratch", # NEED TO MODIFY THESE FOR HPC WORK
    prefix,
    "output",
    "qc_cnvs",
    cnv_qc_sum
  ),
  header = TRUE
)

cnv_include <- read.table(
  file = file.path(
    "./tmp_scratch",
    prefix,
    "output",
    "qc_cnvs",
    cnv_qc_include
  ),
  header = TRUE
)

# cnv_qual=read.table("~/CNV_data/CNV_repo/tmp_scratch/January2026_updated_data_/output/qc_cnvs/January2026_updated_data.qcsum", sep="", header =T)
# cnv_include=read.table("~/CNV_data/CNV_repo/tmp_scratch/January2026_updated_data_/output/qc_cnvs/January2026_updated_data.qcpass", sep="", header =T)


################################################################################
### We get a list of indiviuals who have failed and passed QC based on our predefined parameters
exclude.individuals <- cnv_qual[
  which(as.numeric(as.character(cnv_qual$NumCNV)) >= NCNV_thres |
        as.numeric(as.character(cnv_qual$WF)) >= WF_thres |
        as.numeric(as.character(cnv_qual$WF))<=-WF_thres |
        as.numeric(as.character(cnv_qual$LRR_SD)) >=LRR_SD_thres),]
exclude.individuals$GROUP="QC Fail"

include.individuals <- cnv_qual[
  -which(as.numeric(as.character(cnv_qual$NumCNV)) >= NCNV_thres |
         as.numeric(as.character(cnv_qual$WF)) >=WF_thres |
         as.numeric(as.character(cnv_qual$WF))<=-WF_thres |
         as.numeric(as.character(cnv_qual$LRR_SD)) >=LRR_SD_thres),]
include.individuals$GROUP="QC Pass"

cnv_qual <- rbind(exclude.individuals,include.individuals)

# ==============================================================================
# CNV outliers analysis, visual check
# ==============================================================================
waveform_plot <- ggplot(cnv_qual,aes(x=WF,y=NumCNV,color=GROUP))+
  geom_point(shape=1)+
  xlab("Wavefactor")+
  ylab("N CNVs")+
  scale_color_manual(values=c("#E69F00","#56B4E9"))+
  theme(legend.position="none")

LRR_plot <- ggplot(cnv_qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
  geom_point(shape=1)+
  xlab("LRR SD")+
  ylab("N CNVs")+
  scale_color_manual(values=c("#E69F00","#56B4E9"))+
  theme(legend.position="bottom")


baf_plot <- ggplot(cnv_qual, aes(x=BAF_drift)) + geom_histogram(alpha=.5, position="identity",colour="black", fill="white")+
  xlab("B-Allele Frequency Drift")+
  ylab("Frequency")

LRR_SD_freq_plot <- ggplot(cnv_qual, aes(x=LRR_SD)) + geom_histogram(colour="black", fill="white")+
  xlab("LRR SD")+
  ylab("Frequency")

no_outliers_removed <- grid.arrange(arrangeGrob(waveform_plot,LRR_plot,ncol=1),arrangeGrob(baf_plot,LRR_SD_freq_plot,ncol=1),ncol=2, widths=c(1,1))

# ggsave(plot = no_outliers_removed, 
#        filename = "./tmp_scratch/")

# TODO save this plot as a figure in a figs dir.

################################################################################

LRR_SD_vd <- cnv_qual$LRR_SD >= LRR_SD_thres
NCNV_vd <- cnv_qual$NumCNV >= NCNV_thres
WF_vd <- (cnv_qual$WF <= (-1 * WF_thres) |  cnv_qual$WF >= WF_thres)


vd_all <- cbind(LRR_SD_vd, NCNV_vd, WF_vd)
vdiagram_lowest_thres <- vennCounts(vd_all)

vennDiagram(vd_all, names = c(paste("LRR_SD"), paste("NCNV"), paste("WF")))
title(paste("Individuals excluded using the following parameters: \nLRR>", LRR_SD_thres," ; NCNVs>",NCNV_thres," ; WF>", WF_thres,sep=""))

#######################################################################################

exclude.individuals_ = cnv_qual[which(as.numeric(as.character(cnv_qual$NumCNV)) >= NCNV_thres | as.numeric(as.character(cnv_qual$WF)) >=WF_thres | as.numeric(as.character(cnv_qual$WF))<=-WF_thres | as.numeric(as.character(cnv_qual$LRR_SD)) >=LRR_SD_thres),]
cnv_qual_ = cnv_qual[! cnv_qual$File %in% unlist(exclude.individuals$File),]


# Can probably use this function above PR
plot_cnv_vs_metric <- function(xvar, group_length = 1){

  # map input to pretty label
  xlabels <- c(
    "WF" = "Wave Factor",
    "LRR_SD" = "LRR SD"
  )

  plot=ggplot(cnv_qual,aes(x=.data[[xvar]],y=NumCNV,color=GROUP)) +
    geom_point(shape=1) +
    xlab(xlabels[[xvar]]) +
    ylab("N CNVs")
    if (group_length == 1){
      plot <- plot +
      scale_color_manual(values=c("#56B4E9"))+
        theme(legend.position="none")
    }
    else if (group_length == 2){
      plot <- plot +
      scale_color_manual(values=c("#E69F00", "#56B4E9")) +
        theme(legend.position="none")
    }
  return (plot)
}

if (length(unique(cnv_qual$GROUP)) == 1){
  plot_1 <- plot_cnv_vs_metric("WF", group_length = 1)
}else if(length(unique(cnv_qual$GROUP))==2){
  plot_1 <- plot_cnv_vs_metric("WF", group_length = 2)
}

if (length(unique(cnv_qual$GROUP)) == 1){
  plot_2 <- plot_cnv_vs_metric("LRR_SD", group_length = 1)
}else if(length(unique(cnv_qual$GROUP))==2){
  plot_2 <- plot_cnv_vs_metric("LRR_SD", group_length = 2)
}

# if(length(unique(cnv_qual$GROUP))==1){
#   plot1no=ggplot(cnv_qual,aes(x=WF,y=NumCNV,color=GROUP)) +
#     geom_point(shape=1) +
#     xlab("Wavefactor") +
#     ylab("N CNVs") +
#     scale_color_manual(values=c("#56B4E9")) +
#     theme(legend.position="none")
# } else if(length(unique(cnv_qual$GROUP))==2){
#   plot1no=ggplot(cnv_qual,aes(x=WF,y=NumCNV,color=GROUP)) +
#     geom_point(shape=1) +
#     xlab("Wavefactor") +
#     ylab("N CNVs") +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     theme(legend.position="none")
# }

# if(length(unique(cnv_qual$GROUP))==1){
#   plot2no=ggplot(cnv_qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
#     geom_point(shape=1)+
#     xlab("LRR SD")+
#     ylab("N CNVs")+
#     scale_color_manual(values=c("#56B4E9"))+
#     theme(legend.position="bottom")
# } else if(length(unique(cnv_qual$GROUP))==2){
#   plot2no=ggplot(cnv_qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
#     geom_point(shape=1)+
#     xlab("LRR SD")+
#     ylab("N CNVs")+
#     scale_color_manual(values=c("#E69F00","#56B4E9"))+
#     theme(legend.position="bottom")
# }

# if(length(unique(cnv_qual$GROUP))==1){
#   plot2no=ggplot(cnv_qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
#     geom_point(shape=1)+
#     xlab("LRR SD")+
#     ylab("N CNVs")+
#     scale_color_manual(values=c("#56B4E9"))+
#     theme(legend.position="bottom")
# } else if(length(unique(cnv_qual$GROUP))==2){
#   plot2no=ggplot(cnv_qual,aes(x=LRR_SD,y=NumCNV,color=GROUP))+
#     geom_point(shape=1)+
#     xlab("LRR SD")+
#     ylab("N CNVs")+
#     scale_color_manual(values=c("#E69F00","#56B4E9"))+
#     theme(legend.position="bottom")
# }
# Probably can put this in a fucntion. But only 1 repeat so far. Check rest of script PR

plot_3 <- ggplot(cnv_qual, aes(x=BAF_drift)) + geom_histogram(alpha=.5, position="identity",colour="black", fill="white")+
  xlab("B-Allele Frequency Drift")+
  ylab("Frequency")
#geom_vline(xintercept=sd(cnv_qual[which(cnv_qual$GROUP=="PennCNV Pass"),]$BAF_drift)*2, linetype="dashed", color = "Red")

plot_4 <- ggplot(cnv_qual, aes(x=LRR_SD)) + geom_histogram(colour="black", fill="white")+
  xlab("LRR SD")+
  ylab("Frequency")

grid.arrange(arrangeGrob(plot_1, plot_2,ncol=1),arrangeGrob(plot_3, plot_4,ncol=1),ncol=2, widths=c(1,1))
#############################################

###added by Jo H to combine datafiles
# TODO make this dynamic for name from yaml file!!



# Get a list of all file names in the folder (e.g., CSV files)
file_list <- list.files(path = FOLDER_PATH, pattern = "\\.goodcnv$", full.names = TRUE)

# Use lapply to read all files with fread and bind them together
combined_goodcnv <- rbindlist(lapply(file_list, fread, header = FALSE, fill = TRUE))

###################################################################################################

### Sample /CNV exclusions

CNV_Calls <- as.data.frame(combined_goodcnv)
names(CNV_Calls) <- c("COORDS","NPROBES","SIZE","TYPE","ID","START_PROBE","END_PROBE","CONF")

CNV_Calls$NPROBES=gsub("numsnp=","",CNV_Calls$NPROBES)
CNV_Calls$SIZE=gsub("length=","",CNV_Calls$SIZE)
CNV_Calls$SIZE=gsub(",","",CNV_Calls$SIZE)
CNV_Calls$CHR=gsub("chr","",do.call(rbind,strsplit(CNV_Calls$COORDS,split=":"))[,1])
CNV_Calls$CONF=gsub("conf=","",CNV_Calls$CONF)

CNV_Calls$CONF=as.numeric(as.character(CNV_Calls$CONF)) # Why doing this? They are float numbers? If NA (or something else) we should remove
CNV_Calls$START=as.numeric(gsub(",","",do.call(rbind,strsplit(gsub("chr","",do.call(rbind,strsplit(CNV_Calls$COORDS,split="-"))[,1]),split=":"))[,2]))
CNV_Calls$END=as.numeric(gsub(",","",do.call(rbind,strsplit(CNV_Calls$COORDS,split="-"))[,2]))

CNV_Calls <- CNV_Calls[which(as.numeric(as.character(CNV_Calls$SIZE))>100000) & as.numeric(as.character(CNV_Calls$NPROBES))>=20 & CNV_Calls$CONF>=10,]
#CNV_Calls_Exclude=CNV_Calls[which(as.numeric(as.character(CNV_Calls$SIZE))<20000 | as.numeric(as.character(CNV_Calls$NPROBES))<10 | CNV_Calls$CONF<10),]
CNV_Calls_Exclude <- CNV_Calls[which(as.numeric(as.character(CNV_Calls$SIZE))<100000 | as.numeric(as.character(CNV_Calls$NPROBES))<10 | CNV_Calls$CONF<10),]

#CNV_Calls=CNV_Calls[grep("_R1$",CNV_Calls$ID),]

CNV_Calls_HQ <- CNV_Calls[!CNV_Calls$ID %in% unlist(exclude.individuals$File),] ### this is file with QC fails excluded I think
CNV_Calls_QCFail <- CNV_Calls[CNV_Calls$ID %in% unlist(exclude.individuals$File),]##Jo H

hundredKB_ <- CNV_Calls[which(as.numeric(CNV_Calls_HQ$SIZE) > 100000),] #???
hundredKB <- CNV_Calls # ??? PR 

########neuroCNVS##########################################################################

cnv_neurodev <- read.table(file=cnv_kk_file, sep="\t", stringsAsFactors=F)

cnv_neurodev$CHR <- gsub("chr","",do.call(rbind,strsplit(cnv_neurodev$V3,split=":"))[,1])
cnv_neurodev$START <- as.numeric(gsub(",","",do.call(rbind,strsplit(gsub("chr","",do.call(rbind,strsplit(cnv_neurodev$V3,split="-"))[,1]),split=":"))[,2]))
cnv_neurodev$END <- as.numeric(gsub(",","",do.call(rbind,strsplit(cnv_neurodev$V3,split="-"))[,2]))
cnv_neurodev$CLASS <- do.call(rbind,lapply(1:nrow(cnv_neurodev),function(x){
  ifelse(length(grep("del",cnv_neurodev[x,1]))>0,yes=1,no=3)
}))
names(cnv_neurodev)[2] <- c("LOCUS")
cnv_neurodev <- cnv_neurodev[,c(1,4:7)]

##################################################################################

cnv_neuro_beta <- unique(as.data.frame(do.call(rbind,lapply(1:nrow(CNV_Calls_HQ),function(x){
  cnv_chr <- CNV_Calls_HQ[x,9]
  cnv_start <- CNV_Calls_HQ[x,10]
  cnv_end <- CNV_Calls_HQ[x,11]

  ### CNV spans entire region
  cnv_neurodev.match <- cnv_neurodev[which(cnv_neurodev$CHR==cnv_chr & cnv_start<=cnv_neurodev$START & cnv_end>=cnv_neurodev$END),]

  if(nrow(cnv_neurodev.match)>0){
    cnv.prop=1
    cnv_neurodev.match=cbind(cnv_neurodev.match,cnv.prop)
  } else if(nrow(cnv_neurodev.match)==0){

    ### CNV start is less than the start of the nd cnv and end of the CNV is less than than end of the nd cnv
    cnv_neurodev.match=cnv_neurodev[which(cnv_neurodev$CHR==cnv_chr & cnv_start<=cnv_neurodev$START & cnv_end<=cnv_neurodev$END & cnv_end>=cnv_neurodev$START),]

    if(nrow(cnv_neurodev.match)>0){
      cnv_neurodev.match.start=cnv_neurodev.match$START
      cnv_neurodev.match.end=cnv_neurodev.match$END

      cnv.prop=(cnv_end-cnv_neurodev.match$START)/(cnv_neurodev.match.end-cnv_neurodev.match.start)
      cnv_neurodev.match=cbind(cnv_neurodev.match,cnv.prop)

    }
    ### CNV start is greater than start of nd cnv and ends after the end of the nd region
    else if(nrow(cnv_neurodev.match)==0){
      cnv_neurodev.match=cnv_neurodev[which(cnv_neurodev$CHR==cnv_chr & cnv_start>=cnv_neurodev$START & cnv_end>=cnv_neurodev$END & cnv_start<=cnv_neurodev$END),]

      if(nrow(cnv_neurodev.match)>0){
        cnv_neurodev.match.start=cnv_neurodev.match$START
        cnv_neurodev.match.end=cnv_neurodev.match$END

        cnv.prop=(cnv_neurodev.match$END-cnv_start)/(cnv_neurodev.match.end-cnv_neurodev.match.start)
        cnv_neurodev.match=cbind(cnv_neurodev.match,cnv.prop)

      }
      ### CNV falls inside nd region
      else if(nrow(cnv_neurodev.match)==0){

        cnv_neurodev.match=cnv_neurodev[which(cnv_neurodev$CHR==cnv_chr & cnv_start>=cnv_neurodev$START & cnv_end<=cnv_neurodev$END),]

        if(nrow(cnv_neurodev.match)>0){
          cnv_neurodev.match.start=cnv_neurodev.match$START
          cnv_neurodev.match.end=cnv_neurodev.match$END

          cnv.prop=(cnv_end-cnv_start)/(cnv_neurodev.match.end-cnv_neurodev.match.start)
          cnv_neurodev.match=cbind(cnv_neurodev.match,cnv.prop)
        }}}}

  if(nrow(cnv_neurodev.match)>0){
    cbind(CNV_Calls_HQ[x,],cnv_neurodev.match)
  }

}))))

cnv_neuro_beta$TYPE <- do.call(rbind,strsplit(cnv_neuro_beta$TYPE,split="="))[,2]
cnv_neuro_beta_fifpc <- cnv_neuro_beta[which(cnv_neuro_beta$TYPE==cnv_neuro_beta$CLASS),]

#cnv_neuro_beta.fifpc=cnv_neuro_beta # PR ??????




###################################################################################################

#geneloc=read.table(file=paste("~/Hybrid.Magma/Data/NCBI37.3.gene.loc",sep=""),header=F,stringsAsFactors = F)
geneloc=read.table(file="./exon_files/NCBI37.3.gene.loc",
                     ,header=F,stringsAsFactors = F)


cnv_neuro_beta_genes <- as.data.frame(do.call(rbind,lapply(1:nrow(cnv_neuro_beta_fifpc),function(x){
  a=cnv_neuro_beta_fifpc[x,]
  cnv_chr=a[1,9]
  cnv_start=a[1,10]
  cnv_end=a[1,11]

  geneloc.chr=geneloc[which(geneloc$V2==cnv_chr),]

  ### CNV covers whole gene
  gene.whole=geneloc.chr[which(cnv_start<geneloc.chr$V3 & cnv_end>geneloc.chr$V4),]
  gene.start=geneloc.chr[which(cnv_start<geneloc.chr$V3 & cnv_end<geneloc.chr$V4 & cnv_end>geneloc.chr$V3),]
  gene.end=geneloc.chr[which(cnv_start>geneloc.chr$V3 & cnv_end>geneloc.chr$V4 & cnv_start<geneloc.chr$V4),]
  gene.mid=geneloc.chr[which(cnv_start>geneloc.chr$V3 & cnv_end<geneloc.chr$V4),]


  if(nrow(gene.whole)>0 | nrow(gene.start)>0 | nrow(gene.end)>0 | nrow(gene.mid)>0){
    a$GENES=paste(unique(c(gene.whole$V6,gene.start$V6,gene.end$V6,gene.mid$V6)),collapse="|")
    a
  } else{
    a$GENES=""
    a
  }
})))

##########################################################

Checks<- CNV_Calls_HQ%>%
  dplyr::filter(!ID %in% cnv_neuro_beta_genes$ID) # need to use dplyr:: here?

### Check for duplicate IDs. It would be unusual for an individual to have two different ND CNVs
### Most likely cause for dupliucates are having a smaller nested CNV within a larger one, where both are ND
### Loop over the IDs, check to see if the CNV coordinates are duplicates, and if so just report the largest CNV.

#if(nrow(cnv_neuro_beta_genes[which(duplicated(cnv_neuro_beta_genes$ID)),])>0){

#nested.cnv.dups=as.data.frame(do.call(rbind,lapply(cnv_neuro_beta_genes[which(duplicated(cnv_neuro_beta_genes$ID)),]$ID,function(x){
#  cnv.dups=cnv_neuro_beta_genes[which(cnv_neuro_beta_genes$ID==x),]

#  if(length(unique(cnv.dups$COORDS))==1){

### Get the smaller of the two CNVs (the smaller one will be the nested one)
#    cnv.nested=which.min(cnv.dups[,15]-cnv.dups[,14])
#    cnv.exclude=cnv.dups[cnv.nested,c(5,12)]
#  }
#})))
#}

#if(exists("nested.cnv.dups")==T){
#cnv_neuro_beta.fifpc=cnv_neuro_beta_genes[!(cnv_neuro_beta_genes$ID %in% unlist(nested.cnv.dups[,1]) & cnv_neuro_beta_genes$V1 %in% unlist(nested.cnv.dups[,2])),]
#}

#miscfolder="~/Downloads/"

gene_df <- data.frame(
  process_name = c("nrxn1", "ywhae", "pafah1b1"),
  write_name   = c("NRXN1 Exon", "YWHAE Exon", "PAFAH1B1 Exon "),
  file_address = c("./exon_files/exon.NCBI 1.NRXN1",
                   "./exon_files/exon.YWHAE",
                   "./exon_files/exon.PAFAH1B1"),
  stringsAsFactors = FALSE
)

# Improve address detection here
exon_processing <- function(gene){
  #browser()
  gene_table <- read.table(file=gene_df$file_address[gene_df$process_name == gene], header=T, stringsAsFactors = F)
  gene_exon_start <- unlist(strsplit(gene_table$exonStarts,split = ","))
  gene_exon_end <- unlist(strsplit(gene_table$exonEnds,split = ","))
  gene_exon_coords <- unique(as.data.frame(do.call(rbind,lapply(1:length(gene_exon_start),function(x){
    gene_exon_start_loop = gene_exon_start[x]
    gene_exon_end_loop=gene_exon_end[x]
    cbind(gene_exon_start_loop,gene_exon_end_loop)
  }))))
gene_exon_coords$EXON=paste(gene_df$write_name[gene_df$process_name == gene] ,seq(1:nrow(gene_exon_coords)),sep="")
}


#### NRXN1 exon processing ### NO IDEA WHY THIS ISNT WORKING! TODO check this and get it working. THough not essential
#nrxn1_file <- exon_processing("nrxn1")
#ywhae_file <- exon_processing("ywhae")#  Something is wrong with the function when processing ywhae. Not sure what but can just run manually below for now PR
#pafah_file <- exon_processing("pafah1b1") # So strange, if I comment the line above out, pfah1b1 fails....

nrxn1=read.table(file="./exon_files/exon.NCBI1.NRXN1",header=T,stringsAsFactors = F)
nrx1.exon.start=unlist(strsplit(nrxn1$exonStarts,split = ","))
nrx1.exon.end=unlist(strsplit(nrxn1$exonEnds,split = ","))
nrx1.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(nrx1.exon.start),function(x){
  nrx1.exon.start.loop=nrx1.exon.start[x]
  nrx1.exon.end.loop=nrx1.exon.end[x]
  cbind(nrx1.exon.start.loop,nrx1.exon.end.loop)
}))))
nrx1.exon.coords$EXON=paste("NRXN1 Exon ",seq(1:nrow(nrx1.exon.coords)),sep="")
###
#
### YWHAE exon processing #17p
YWHAE=read.table(file="./exon_files/exon.YWHAE",header=T,stringsAsFactors = F,fill = T)
YWHAE.exon.start=unlist(strsplit(YWHAE$exonStarts,split = ","))
YWHAE.exon.end=unlist(strsplit(YWHAE$exonEnds,split = ","))
YWHAE.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(YWHAE.exon.start),function(x){
  YWHAE.exon.start.loop=YWHAE.exon.start[x]
  YWHAE.exon.end.loop=YWHAE.exon.end[x]
  cbind(YWHAE.exon.start.loop,YWHAE.exon.end.loop)
}))))
YWHAE.exon.coords$EXON=paste("YWHAE Exon ",seq(1:nrow(YWHAE.exon.coords)),sep="")
###
#
### PAFAH1B1 exon processing #17p
PAFAH1B1=read.table(file="./exon_files/exon.PAFAH1B1",header=T,stringsAsFactors = F,fill = T)
PAFAH1B1.exon.start=unlist(strsplit(PAFAH1B1$exonStarts,split = ","))
PAFAH1B1.exon.end=unlist(strsplit(PAFAH1B1$exonEnds,split = ","))
PAFAH1B1.exon.coords=unique(as.data.frame(do.call(rbind,lapply(1:length(PAFAH1B1.exon.start),function(x){
  PAFAH1B1.exon.start.loop=PAFAH1B1.exon.start[x]
  PAFAH1B1.exon.end.loop=PAFAH1B1.exon.end[x]
  cbind(PAFAH1B1.exon.start.loop,PAFAH1B1.exon.end.loop)
}))))
PAFAH1B1.exon.coords$EXON=paste("PAFAH1B1 Exon ",seq(1:nrow(PAFAH1B1.exon.coords)),sep="")
###

cnv_neuro_beta_genes$CRITERIA_MET=0
cnv_neuro_beta_genes$SIZE=as.numeric(cnv_neuro_beta_genes$SIZE)

cnvs_unique <- unique(cnv_neuro_beta_genes$V1)


# This is assigning the criteria met column to 1 or 0 depending on a number of conditions using if else statements
cnv_patho_criteria <- as.data.frame(do.call(rbind,lapply(cnvs_unique,function(x){
  # browser()
  #cnv_patho_criteria[1,12]



  CNV.ND=cnv_neuro_beta_genes[which(cnv_neuro_beta_genes$V1==x),]

  as.data.frame(do.call(rbind,lapply(1:nrow(CNV.ND),function(y){
    CNV.ND.Line=CNV.ND[y,]

    greater_than_1mbp <- CNV.ND.Line$SIZE > 1000000
    greater_than_4mbp <- CNV.ND.Line$SIZE > 4000000
    size_greater_50 <- CNV.ND.Line$cnv.prop > 0.5
    size_greater_80 <- CNV.ND.Line$cnv.prop > 0.8

    # 1q21.1	del/dup		 Size	>50%	of	critical	region
    if((x=="1q21.1 del" | x=="1q21.1 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 1p36 del/dup		 Size	>50%	of	critical	region,	affecting	GABRD
    if((x=="1p36 del (GABRD)" | x=="1p36 dup (GABRD)") & length(grep("GABRD",CNV.ND.Line$GENES))>0 & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # TAR del/dup		 Size	>50%	of	critical	region
    if((x=="TAR del" | x=="TAR dup") & size_greater_50){
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
    if((x=="2q11.2 del" | x=="2q11.2 dup") & length(grep("LMAN2L",CNV.ND.Line$GENES))>0 & length(grep("ARID5A",CNV.ND.Line$GENES))>0 & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 2q13	del/dup		 Size	>50%	of	critical	region
    if((x=="2q13 del" | x=="2q13 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 2q13	del/dup	(NPHP1)		 Size	>50%	of	critical	region,	affecting	NPHP1
    if((x=="2q13 del (NPHP1)" | x=="2q13 dup (NPHP1)") & length(grep("NPHP1",CNV.ND.Line$GENES))>0 & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 2q21.1	del/dup		 Size	>50%	of	critical	region
    if((x=="2q21.1 del" | x=="2q21.1 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 2q37	del	(HDAC4)		 Size	>50%	of	critical	region,	affecting	HDAC4
    if((x=="2q37 del (HDAC4)" | x=="2q37 dup (HDAC4)") & size_greater_50 & length(grep("HDAC4",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 3q29	del/dup		 Size	>50%	of	critical	region
    if((x=="3q29	del" | x=="3q29	dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # WolfâHirschhorn	del/dup		 Size	>50%	of	critical	region
    if((x=="WolfâHirschhor del" | x=="WolfâHirschhor dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # Sotos	Syn/5q35	dup		 Size	>50%	of	critical	region
    if(x=="Sotos Syn/5q35 dup" & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # Williams	Beuren	Syn	del/dup	Size	>50%	of	critical	region
    if((x=="Williams	Beuren Syn del" | x=="Williams Beuren Syn dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 8p23.1	del/dup		 At	least	1Mbp	of	critical	region
    if((x=="8p23.1 del" | x=="8p23.1 dup") & greater_than_1mbp){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 9q34	del/dup	(EHMT1)		 At	least	1Mbp	CNVs,	including	EHMT1
    if(x=="9q34 del (EHMT1)" & greater_than_1mbp & length(grep("EHMT1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 10q11.21q11.23	del/dup		 Size	>50%	of	critical	region
    if((x=="10q11.21q11.23 del" | x=="10q11.21q11.23 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 10q23	del/dup		 At	least	1Mbp,	including	NRG3 and	GRID1
    if((x=="10q23 del" | x=="10q23 dup") & greater_than_1mbp & length(grep("NRG3",CNV.ND.Line$GENES))>0 & length(grep("GRID1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # PotockiâShaffer	Syn	del/11p11.2	dup (EXT2), size	>50%	of	critical	region,	including	EXT2
    if(x=="Potocki-Shaffer syndrome del (EXT2)" & size_greater_50 & length(grep("EXT2",CNV.ND.Line$GENES))){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 15q11.2	del/dup		 Size	>50%	of	critical	region
    if((x=="15q11.2 del" | x=="15q11.2 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 15q11.2 del/dup BP1-BP2		 Size	>50%	of	critical	region
    if((x=="15q11.2 del BP1-BP2" | x=="15q11.2 dup BP1-BP2") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # PWS	del/dup		 Full	critical	region,	~4Mbp
    if((x=="(PWS/AS del" | x=="(PWS/AS dup") & greater_than_4mbp & CNV.ND.Line$cnv.prop>0.8){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 15q13.3 del BP4-BP5		 Size	>50%	of	critical	region
    if(x=="15q13.3 del BP4-BP5" & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 15q24	del/dup		 At	least	1Mbp	between	the	AâE	intervals
    if((x=="15q24 del" | x=="15q24 dup") & greater_than_1mbp){
      CNV.ND.Line$CRITERIA_MET=1
    }



    # 15q25	del/dup		 At	least	1Mbp	between	the	AâD	intervals
    if(x=="15q25 del" & greater_than_1mbp){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 16p13.11	del/dup		 Size	>50%	of	critical	region
    if((x=="16p13.11 del" | x=="16p13.11 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 16p12.1	del		 Size	>50%	of	critical	region
    if(x=="16p12.1 del" & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }



    # 16p11.2	distal	del/distal	dup		 Size	>50%	of	critical	region
    if((x=="16p11.2 distal del" | x=="16p11.2 distal dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 16p11.2	del/dup		 Size	>50%	of	critical	region
    if((x=="16p11.2 del" | x=="16p11.2 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }



    # 17p13.3	del/dup	(YWHAE)		 Exonic	deletions;	whole	gene	duplications

    if((x=="17p13.3 del (YWHAE)" | x=="17p13.3	dup (YWHAE") & size_greater_50 & length(grep("YWHAE",CNV.ND.Line$GENES))>0){

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
    if((x=="17p13.3 del (PAFAH1B1)" | x=="17p13.3	dup (PAFAH1B1") & size_greater_50 & length(grep("PAFAH1B1",CNV.ND.Line$GENES))>0){
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
    if((x=="17q11.2 del (NF1)" | x=="17q11.2 dup (NF1)") & size_greater_50 & length(grep("NF1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1
    }


    # SmithâMagenis/PotockiâLupski	Syndrome	Size	>50%	of	critical	region

    if((x=="Smith-Magenis syndrome del" | x=="Potocki-Lupski syndrome dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 17q11.2	del/dup	(NF1)		 Size	>50%	of	critical	region,	affecting	NF1
    if((x=="17q11.2 del (NF1)" | x=="17q11.2 dup (NF1)") & size_greater_50 & length(grep("NF1",CNV.ND.Line$GENES))>0){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 17q12	del/dup		 Size	>50%	of	critical	region
    if((x=="17q12 del" | x=="17q12 dup" | x=="Renal cysts and diabetes syndrome del") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 17q21.31	del/dup		 Size	>50%	of	critical	region
    if((x=="17q21.31 del" | x=="17q21.31 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 17q23.1q23.2	del		 Size	>50%	of	critical	region
    if(x=="17q23.1q23.2 del" & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 17p12       del              Size   >50%    of      critical        region
    if(x=="17p12 del" & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 22q11.2	del/dup		 Size	>50%	of	critical	region
    if((x=="22q11.2 del" | x=="22q11.2 dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # 22q11.2	distal	del/dup		 Size	>50%	of	critical	region
    if((x=="22q11.2 distal del" | x=="22q11.2 distal dup") & size_greater_50){
      CNV.ND.Line$CRITERIA_MET=1
    }

    # SHANK3 del/dup	 At	least	1Mbp	CNVs,	including	SHANK3
    if((x=="SHANK3 del" | x=="SHANK3 dup") & length(grep("SHANK3",CNV.ND.Line$GENES))>0 & greater_than_1mbp){
      CNV.ND.Line$CRITERIA_MET=1
    }


    CNV.ND.Line

  })))
})))

patho_criteria_met <- cnv_patho_criteria[which(cnv_patho_criteria$CRITERIA_MET==1),]

select_non_nested_cnvs <- function(x, dataset){
  
  #browser()
  
  # Select rows for this individual
  a <- dataset[dataset$ID == x, ]
  
  #Single row, nothing to filter
  if (nrow(a) == 1){
    return(a)
  }
  
  # Multiple rows, initialise empty outputs 
  b1 <- b2 <- b3 <- NULL
  
  # Rule1 : TAR nested in 1q21
  if ((length(grep("TAR", a$V1)) + length(grep("1q21", a$V1))) == 2) {
    b1 <- a[-grep("TAR", a$V1), ]
  }
  
  # Rule2 : 15q11.2 nested in PWS/AS
  if ((length(grep("15q", a$V1)) + length(grep("PWS", a$V1))) == 2) {
    b2 <- a[-grep("15q", a$V1), ]
  }
  
  # Rule3 : 16p11.2 nested 
  if ((length(grep("distal", a$V1)) + length(grep("16p11.2 del", a$V1))) == 2) {
    b3 <- a[-grep("distal", a$V1), ]
  }
  
  # Final fallback 
  
  if (!is.null(b1)) return(b1)
  if (!is.null(b2)) return(b2)
  if (!is.null(b3)) return(b3)
  
  # Nothing matched → return all rows
  return(a)
}



patho_criteria_met.no.nested <- as.data.frame(
  do.call(
    rbind,
    lapply(
      unique(patho_criteria_met$ID),
      select_non_nested_cnvs,
      dataset = patho_criteria_met
    )
  )
)
### Check I do this right. I still dont really understand this PR 
write.csv(patho_criteria_met, file = "./tmp_scratch/January2026_updated_data/output/Routput/patho_criteria_met.csv")
write.csv(patho_criteria_met.no.nested, file = "./tmp_scratch/January2026_updated_data/output/Routput/patho_criteria_met_no_nested.csv")

patho_criteria_met = patho_criteria_met.no.nested # Why? PR

patho_criteria_met = patho_criteria_met[order(as.numeric(as.character(patho_criteria_met[,9])),as.numeric(as.character(patho_criteria_met[,10]))),]
cnv.patho=as.data.frame(table(patho_criteria_met$V1))
cnv.patho=unique(merge(patho_criteria_met[,c(12,9,10)],cnv.patho,by.x="V1",by.y="Var1",sort=F)[,c(1,4)])


cnv.patho$DATASET="JAN26"
names(cnv.patho)=c("Neurodevelopmental CNV","N","DATASET")

CalledCNVS<- rbind(patho_criteria_met.no.nested,cnv_patho_criteria[which(cnv_patho_criteria$CRITERIA_MET==0),])
# TODO for now this should work, but need to add this into the YAML to make seamless.

parent_address <- paste0("./tmp_scratch/", prefix, "/output/") # TODO change this
subdir_routput <- "Routput"
dir.create(file.path(parent_address, subdir_routput), showWarnings = FALSE)



write.table(CalledCNVS,
            file = paste0("./tmp_scratch/", prefix, "/output/Routput/called_cnvs.txt"),
            col.names=T,
            row.names=F,
            quote=F,
            sep="\t")
write.csv(CalledCNVS,
            file = paste0("./tmp_scratch/", prefix, "/output/Routput/called_cnvs.csv"),
            col.names=T)

# ==============================================================================
# Generate the split plots directory 
# ==============================================================================
split_file_dir <- paste0("./tmp_scratch/", prefix, "/output/split_files/")

r_output_address <- paste0("./tmp_scratch/", prefix, "/output/Routput/")
subdir_plots <- "plots"
dir.create(file.path(r_output_address, subdir_plots), showWarnings = FALSE)


# A lot of this is copy paste. Make into a simple check
# TODO get this working, ideally would then add a seperate text file with all required information. Perhaps still the HTML file if needed.
generate_CNV_plots <- function(dataset){

}


# cnv.raw.params

invisible(lapply(unique(patho_criteria_met$ID),function(x){
  #browser()
  # This is taking the whole row of unique ID
  patho.id=patho_criteria_met[which(patho_criteria_met$ID==x),]
  
  magic_number_1 <- 1.3 ## Ask what this means?
  # Helper function to grab coordinates, used in both if and else statements below. We should ask about the 1.3 magic number and what this does? 
  get_coords <- function(patho_id_file){

    patho.cnv.lower.coords <- as.numeric(as.character(patho_id_file[,10]))-(patho_id_file[1,3]* magic_number_1)
    patho.cnv.upper.coords <- as.numeric(as.character(patho_id_file[,11]))+(patho_id_file[1,3]* magic_number_1)
    patho_cnv_chr <- as.numeric(as.character(patho_id_file[,9]))
    
    coords_list <- list(
      lower_coords = patho.cnv.lower.coords,
      upper_coords = patho.cnv.upper.coords,
      patho_cnv_chr = patho_cnv_chr
    )
    
    return (coords_list)
  }
  
  if(nrow(patho.id)==1){
    IDname_long <- basename(patho.id$ID)
    IDname <- sub(".*\\.", "", IDname_long)
    coords_list <- get_coords(patho.id)

    cnv.raw <- as.data.frame(fread(input=paste0(parent_address, "/split_files/", x)),header=T,sep="\t")# Reads the filename into fread.  
    
    
    cnv.raw.params <- cnv.raw[
      which(
        cnv.raw$Position >= coords_list[["lower_coords"]] &
          cnv.raw$Position <= coords_list[["upper_coords"]] &
          cnv.raw$Chr == coords_list[["patho_cnv_chr"]]),]
    
    cnv.raw.params$GROUP <- c("Probes outside called CNV")
    cnv.raw.params[
                    which(cnv.raw.params$Position>=as.numeric(as.character(patho.id[,10])) &
                          cnv.raw.params$Position<=as.numeric(as.character(patho.id[,11])) &
                          cnv.raw.params$Chr == coords_list[["patho_cnv_chr"]]),]$GROUP = c("Probes within individually called CNV")
    cnv.raw.params$GROUP <- factor(cnv.raw.params$GROUP,
                                   levels = c("Probes within individually called CNV","Probes outside called CNV"))
    
    baf.col=grep("Allele",names(cnv.raw.params))
    logr.col=grep("Log",names(cnv.raw.params))


    # Generates figures that show BAF and LRR for the split files plots directory

    base_position_figure <- function(data, baf.col, log.col, IDname, type = "baf"){
      #browser()
      if (type == "baf"){
        y_params = data[,baf.col]
        y_title <- "B-Allele Frequency"
        intercept_val <- 0.5
        y_lims <- c(0,1)
      }
      else if (type == "lrr"){
        y_params = data[,logr.col]
        y_title <- "Log R Ratio"
        intercept_val <- 0
        y_lims <- c(-1,1)
      }
      else {
        stop("Unknown type. Use 'baf' or 'lrr'.")
      }
      
      figure=ggplot(data=data,aes(x=data[,3],y=y_params,color=as.factor(GROUP))) +
        geom_point(shape=1) +
        xlab("Base Position")+
        ylab(y_title)+
        ggtitle(paste(IDname, patho.id$V1, sep=" "))+
        geom_hline(yintercept=0.5, linetype="dashed", color = "Black")+
        geom_vline(xintercept=patho.id[,14], linetype="dashed", color = "Green")+
        geom_vline(xintercept=patho.id[,15], linetype="dashed", color = "Green")+
        ylim(y_lims)+
        labs(color="CNV Probe Legend")
      
      return (figure)
    }
    
    baf <- base_position_figure(data = cnv.raw.params, baf.col, baf.log, IDname, "baf")
    lrr <- base_position_figure(data = cnv.raw.params, baf.col, baf.log, IDname, "lrr")

    combined_plot <- grid.arrange(arrangeGrob(baf,lrr,ncol=1))
    
    png(paste0("./tmp_scratch/", prefix, "/output/Routput/plots/",IDname,"___",patho.id$V1,".png",sep=""), width = 10, height = 4, units = 'in', res = 300)
    grid.arrange(arrangeGrob(baf,lrr,ncol=1))
    dev.off()
  }

  else if(nrow(patho.id)>1){
    lapply(1:nrow(patho.id),function(y){
      
      IDname_long <- basename(patho.id$ID)
      IDname <- sub(".*\\.", "", IDname_long)
      coords_list <- get_coords(patho.id)
      
      magic_number_1 <- 1.3 
      

      # Use only the one row for interation
      coords_list <- get_coords(patho.id[y, , drop = FALSE])
      
      cnv.raw <- as.data.frame(fread(input=paste0(parent_address, "/split_files/", x)),header=T,sep="\t")# Reads the filename into fread.  

      cnv.raw.params <- cnv.raw[
        which(
          cnv.raw$Position >= coords_list[["lower_coords"]] &
          cnv.raw$Position <= coords_list[["upper_coords"]] &
          cnv.raw$Chr == coords_list[["patho_cnv_chr"]]),]
      
      cnv.raw.params$GROUP <- c("Probes outside called CNV")
      cnv.raw.params[which(cnv.raw.params$Position>=as.numeric(as.character(patho.id[,10])) & cnv.raw.params$Position<=as.numeric(as.character(patho.id[,11])) & cnv.raw.params$Chr==coords_list[["patho_cnv_chr"]]),]$GROUP=c("Probes within individually called CNV")
      cnv.raw.params$GROUP <- factor(cnv.raw.params$GROUP, levels = c("Probes within individually called CNV","Probes outside called CNV"))

      baf.col=grep("Allele",names(cnv.raw.params))
      logr.col=grep("Log",names(cnv.raw.params))
      
      base_position_figure <- function(data, baf.col, log.col, IDname, type = "baf"){
        #browser()
        if (type == "baf"){
          y_params = data[,baf.col]
          y_title <- "B-Allele Frequency"
          intercept_val <- 0.5
          y_lims <- c(0,1)
        }
        else if (type == "lrr"){
          y_params = data[,logr.col]
          y_title <- "Log R Ratio"
          intercept_val <- 0
          y_lims <- c(-1,1)
        }
        else {
          stop("Unknown type. Use 'baf' or 'lrr'.")
        }
        
        figure=ggplot(data=data,aes(x=data[,3],y=y_params,color=as.factor(GROUP))) +
          geom_point(shape=1) +
          xlab("Base Position")+
          ylab(y_title)+
          ggtitle(paste(IDname, patho.id$V1, sep=" "))+
          geom_hline(yintercept=0.5, linetype="dashed", color = "Black")+
          geom_vline(xintercept=patho.id[,14], linetype="dashed", color = "Green")+
          geom_vline(xintercept=patho.id[,15], linetype="dashed", color = "Green")+
          ylim(y_lims)+
          labs(color="CNV Probe Legend")
        
        return (figure)
      }
      
      baf <- base_position_figure(data = cnv.raw.params, baf.col, baf.log, IDname, "baf")
      lrr <- base_position_figure(data = cnv.raw.params, baf.col, baf.log, IDname, "lrr")
      browser()
      id <- patho.id$V1[1]
      IDname <- IDname[1] # this is just to get a single line for the address path
      combined_plot <- grid.arrange(arrangeGrob(baf,lrr,ncol=1))
      file_path <- paste0(
        "./tmp_scratch/", prefix, "/output/Routput/plots/",
        IDname, "___", id, ".png" # Check with Jo that this is okay? Any edge cases?
      )
      
      if (file.exists(file_path)) {
        file_path <- paste0(
          "./tmp_scratch/", prefix, "/output/Routput/plots/",
          IDname, "___", id, "another_copy.png"
        )
      }
      

      print(file_path)
      
      png(file_path, width = 10, height = 4, units = 'in', res = 300)
      grid.arrange(arrangeGrob(baf,lrr,ncol=1))
      dev.off()
      
    })
  }

}))





# Lets have a go at generating the figure Jess suggests with each breakpoint shown
# We have a LOT of ND-CNVS here so maybe we can make some seperate subsets to make it a bit more readable?

## get out the breakpoints from file

breakpoint_file <- read.csv("./penn_cnv_files/CNV_master_list.txt", header = FALSE)

color_all_cnvs <- ggplot ()


