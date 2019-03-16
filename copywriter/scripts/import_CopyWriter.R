#!/usr/bin/Rscript --vanilla 
# -------------------------------------------------------------------
# Converts CopywriteR output into tables displayable by CNspector
# -------------------------------------------------------------------
library(matrixStats)

input_dir <- "output"
output_dir <- "CNspector"
if(!dir.exists(output_dir)) { dir.create(output_dir,recursive=T) }

display_resolution <- c("targeted","wg_fine","wg_medium","wg_coarse")
bin_resolution <- c("hg19_10kb","hg19_50kb","hg19_500kb","hg19_1000kb")
index_file <- file.path(output_dir,"index.tsv")
counts_cutoff <- 100
df_index <- data.frame(path=c("wg_transcripts.bed","cytoBand.txt","chromInfo.txt","targeted_transcripts.bed"),	
                       file="reference",	
                       sample=c("wg_annotation_file","cytobands_file","chrom_info_file","targeted_annotation_file"),	
                       resolution=c("wg","wg","wg","targeted"),	
                       is_targeted=c(F,F,F,T),	
                       is_wgs=c(T,T,T,F))

names(display_resolution) <- bin_resolution

for (resolution in bin_resolution) {
  df_igv <- read.table(file.path(input_dir,resolution,"CNAprofiles/log2_read_counts.igv"),stringsAsFactors=F,header=T)
  df_read_counts <- read.table(file.path(input_dir,resolution,"CNAprofiles/read_counts.txt"),stringsAsFactors=F,header=T)
  rs <- rowMedians(as.matrix(df_read_counts[,grepl("compensated",colnames(df_read_counts))]),na.rm = T)
  idx_rows <- (rs >= counts_cutoff)
  igv_cols <- colnames(df_igv)
  idx_cols <- grepl("^log2",igv_cols)
  sample_names <- sub("^log2.","",igv_cols[idx_cols])
  colnames(df_read_counts) <- sub(".bam$","",colnames(df_read_counts))
  colnames(df_igv) <- sub(".bam$","",colnames(df_igv))
  sample_names <-  sub(".bam$","",sub("^log2.","",sample_names))
  for (sample in sample_names) {
    df_sample <- cbind(df_read_counts[,c("Chromosome","Start","End")],
                       df_read_counts[,paste0("read.counts.",sample)],
                       df_read_counts[,paste0("read.counts.compensated.",sample)])
    colnames(df_sample) <- c("chr","start","end","raw","cor")
    df_sample$chr <- paste0("chr",df_sample$chr)
    df_sample$N <- NA
    df_sample$N[df_read_counts$Feature %in% df_igv$Feature] <- 2^(1+df_igv[,paste0("log2.",sample)])
    df_sample$Ns <- round(df_sample$N)
    sample_filename <- paste0(sample,"_",resolution,".tsv")
    # Add to index df and write sample file for this resolution
    write.table(df_sample[idx_rows,],file.path(output_dir,sample_filename),sep='\t',quote=F,row.names = F)
    df_index_entry <- data.frame(path=sample_filename,	
                           file="counts",	
                           sample=sample,	
                           resolution=display_resolution[resolution],	
                           is_targeted=(display_resolution[resolution]=="targeted"),	
                           is_wgs=(display_resolution[resolution]!="targeted"))
    df_index <- rbind(df_index,df_index_entry)
  }
}
write.table(df_index,index_file,sep='\t',quote=F,row.names = F)
