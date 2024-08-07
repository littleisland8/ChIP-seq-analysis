.libPaths(new="/home/simone.romagnoli/R/x86_64-pc-linux-gnu-library/4.3/")
library(ChIPQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

rm(list=ls())

## Load sample data
samples <- read.csv('/projects2/2024_Chiarugi/ChIP-seq-analysis/resources/sample_sheet.chipqc.csv')
#View(samples)

blacklist <- read.table("/projects2/2024_Chiarugi/ChIP-seq-analysis/resources/hg38-blacklist.v2.bed", sep="\t", header=TRUE)
blacklist <- makeGRangesFromDataFrame(blacklist,keep.extra.columns = T)

##########  subset to LA vs 10Input
samples_ <- subset(samples,(Factor=="LA_10"))

## Create ChIPQC object
chipObj <- ChIPQC(samples_, annotation="hg38")#, blacklist = "/projects2/2024_Chiarugi/ChIP-seq-analysis/resources/hg38-blacklist.v2.bed") 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: IPLA vs 10_input", reportFolder="/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_10Input")

qc_report <- as.data.frame(QCmetrics(chipObj))
write.csv(qc_report,"/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_10Input/qc_report.csv")


#plotRap(chipObj, facetBy="Condition")


##########  subset to LA vs IgG
samples_ <- subset(samples,(Factor=="LA_IgG"))

## Create ChIPQC object
chipObj <- ChIPQC(samples_, annotation="hg38")#, blacklist = blacklist) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: IPLA vs IgG", reportFolder="/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_IgG")

qc_report <- as.data.frame(QCmetrics(chipObj))
write.csv(qc_report,"/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_IgG/qc_report.csv")

#plotRap(chipObj, facetBy="Condition")


##########  subset to Ctr vs 10Input
samples_ <- subset(samples,(Factor=="CTR_10"))

## Create ChIPQC object
chipObj <- ChIPQC(samples_, annotation="hg38")#, blacklist = blacklist) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: IPCTR vs 10_input", reportFolder="/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPCtr_10Input")

qc_report <- as.data.frame(QCmetrics(chipObj))
write.csv(qc_report,"/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPCtr_10Input/qc_report.csv")

#plotRap(chipObj, facetBy="Condition")


##########  subset to LA vs IgG
samples_ <- subset(samples,(Factor=="CTR_IgG"))

## Create ChIPQC object
chipObj <- ChIPQC(samples_, annotation="hg38", chromosomes = c(paste0("chr",seq(1,22)),"chrX","chrY"))#, blacklist = blacklist) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: IPCTR vs IgG", reportFolder="/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPCtr_IgG")

qc_report <- as.data.frame(QCmetrics(chipObj))
write.csv(qc_report,"/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPCtr_IgG/qc_report.csv")

#plotRap(chipObj, facetBy="Condition")


##########  subset to LA and Ctr vs 10Input
samples <- read.table('/projects2/2024_Chiarugi/ChIP-seq-analysis/resources/sample_sheet.chipqc.tsv', header=TRUE, sep="\t")
samples_ <- subset(samples,(Condition=="LA_10" | Condition=="CTR_10"))

## Create ChIPQC object
chipObj <- ChIPQC(samples_, annotation="hg38", chromosomes = c(paste0("chr",seq(1,22)),"chrX","chrY"))#, blacklist = blacklist) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP QC report: IPLA and IPCtr vs 10 Input", reportFolder="/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_IPCtr_10Input_2")

qc_report <- as.data.frame(QCmetrics(chipObj))
write.csv(qc_report,"/projects2/2024_Chiarugi/ChIP-seq-analysis/qc/ChIPQCreport_IPLA_IPCtr_10Input_2/qc_report.csv")

#plotRap(chipObj, facetBy="Condition")
