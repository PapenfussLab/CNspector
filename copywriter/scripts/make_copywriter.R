#!/usr/bin/Rscript --vanilla
#---------------------------------------
#To run copywriter
#Tested on R 3.5.1 and copywriter 2.10.0
#--------------------------------------

args = commandArgs(trailing=TRUE)

if (length(args)==0) {
	stop("Require the input control samples list", call.=FALSE)}

library(CopywriteR)

#example - copywriter/input/copywriter_sample_controls.txt
input_file<-args[1]
output<-"copywriter_output"
reference<-"reference"
annotations<-"../annotations/browser"

dir.create(reference)
dir.create(output)

bp.param <- SnowParam(workers = 24, type = "FORK", bpprogressbar=T,tasks=24,stop.on.error=F)
bp.param <- MulticoreParam(workers = min(multicoreWorkers(),24),
                           type = "FORK", 
                           bpprogressbar=T,
                           tasks=24,
                           stop.on.error=F)
                           
for (bin_size in c(1000,500,50,10)) {
  sample.control <- read.table(input_file,header=T)
  destination.folder <- paste0(output,"/hg19_",bin_size,"kb")
  reference.folder <- paste0(reference,"/hg19_",bin_size,"kb")
  dir.create(destination.folder)
  capture.regions.file <- paste0(annotations,"/panel_transcripts.bed")
  keep.intermediary.files <- T
  preCopywriteR("reference/", bin_size*1000,"hg19")
  CopywriteR(sample.control=sample.control, 
             destination.folder=destination.folder, 
             reference.folder=reference.folder, 
             bp.param=bp.param,
             capture.regions.file=capture.regions.file, 
             keep.intermediary.files = keep.intermediary.files)  
}
