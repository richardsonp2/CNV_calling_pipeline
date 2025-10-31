# CNV Calling Pipeline

An integrated pipeline for detecting and analyzing copy number variations (CNVs) from SNP array data using Genome Studio, PennCNV, and custom R analysis tools.

## Overview

Much of this file is adapted from the SOP compiled by Alexandra Evans and combines scripts developed by Ellis Pires.

## Prerequisits

I assume you have knowledge of bash, R, a simple overview of cloning and working with git repositories and how to manage filezilla for transfers. The genome studio part of this requires the windows OS, for linux and mac virtual box is a very good package to use. Details about virtualbox installation can be found here.

If you are unfamiliar with bash and shell scripting these resources will be highly useful:

If you are unfamiliar with R these resources are very useful:
<https://r4ds.hadley.nz/>
<https://adv-r.hadley.nz/>
. I code everything I can in R with tidyverse, I try to stay away from base R where I can.

If you need help with filezilla this documentation is useful: <https://wiki.filezilla-project.org/Using>

git gets very complex and complicated very quickly. Gitlab is almost the same as github and the two are interchangable. Just be aware how to clone and how to version control scritpts. A very detailed overview is presented here: <https://git-scm.com/book/en/v2>

## Genome studio

Installation of genome studio is done by following this link:
<https://support.illumina.com/array/array_software/genomestudio/downloads.html>
Choose GenomeStudio v2.0.5 and complete the installation wizard.

Open genome studio in the virtualbox system (or windows if using natively). Genome studio version 2.X or higher (if there is any newer versions...) is required for opening files / projects made in genome studio v2.X. The system is not backwards compatible.

Genome studio recognises .bsc files as Genome studio project files.

When complete
The Genome Studio file is not directly compatible with PennCNV and first needs to be converted into a text file containing only the appropriate columns for each sample.  PennCNV requires the LogR Ratio (LRR) and B Allele Frequency (BAF) as measures of the Signal Intensity data to identify CNVs.

    i. In the Full Data Table, click on the Column Chooser icon
    ![alt text](image.png)

This allows us to choose the columns that are displayed on the Genome Studio file.  PennCNV assumes each sample will contain information in three columns, the genotype, the Log R Ratio and B Allele Frequency.  

    ii. In the Displayed Columns menu, highlight and hide the Index, Address, Gen Train Score, and Frac A/C/T/and G columns. Ensure the Displayed columns only show Name, Chr, Position, and the sample locations.
    iii. Log R Ratio and B Allele Freq will appear in the Hidden Subcolumns window.  Highlight LRR and BAF in this menu and click show to transfer them to the Displayed Subcolumns window.  Ensure the Displayed Subcolumns only shows GType, B Allele Freq and Log R Ratio.

Once this is done, the menu should look like the image below:

Virtual box has a shared folder system. Devices > shared folders > shared folders settings. Mount a folder to a directory on the host. Drag files in the guest and the files appear in the host. Simples🦦.

## Using Penn CNV

Penn CNV is used to **. A handy wrapper is coded in: ```penncnv-calling.sh```
as part of the dpmcn-codebank/shared-pipelines/cnv-calling-using-penncnv/ repository in gitlab.

The command to run uses 4 arguments
the first is the genome studio text file.
The second is the pfb file taken from /nfs/neurocluster/databank/CNV_resources/
The third is the gccode file taken from /nfs/neurocluster/databank/CNV_resources/
the fourth is the main CNV file

```bash
sbatch penncnv-calling.sh /scratch/c.c1837975/cnv_calling/<Genome_Studio_output>.txt /scratch/c.c1837975/cnv_calling/PsychChip_PFB.pfb 
```
