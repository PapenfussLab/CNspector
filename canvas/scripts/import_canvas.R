#!/usr/bin/Rscript --vanilla 
# -------------------------------------------------------------------
# Converts canvas output into tables displayable by CNspector
# -------------------------------------------------------------------
library(genomation)
library(GenomicRanges)
library(stats)
library(VariantAnnotation)

annotations_dir <- "../annotations/preprocessing"
input_dir <- "canvas_output"
output_dir <- "CNspector"
IMPORT_CNV_VCF <- FALSE # The final canvas vcf contains an incomplete set of CNVs and can cause problems.
if(!dir.exists(output_dir)) { dir.create(output_dir,recursive=T) }

index_file <- file.path(output_dir,"index.tsv")
df_index <- data.frame(path=c("wg_transcripts.bed","cytoBand.txt","chromInfo.txt","targeted_transcripts.bed"),	
                       file="reference",	
                       sample=c("wg_annotation_file","cytobands_file","chrom_info_file","targeted_annotation_file"),	
                       resolution=c("wg","wg","wg","targeted"),	
                       is_targeted=c(F,F,F,T),	
                       is_wgs=c(T,T,T,F))
sample_names <- unique(gsub("\\..*","",dir(path=input_dir,pattern="gz$")))
print(sample_names)
for (sample in sample_names) {
  print(sample)
  df_cleaned <- read.table(file.path(input_dir,paste0(sample,".cleaned.gz")),stringsAsFactors=F,header=F,sep='\t')
  df_partitioned_raw <- read.table(file.path(input_dir,paste0(sample,".partitioned.raw.gz")),stringsAsFactors=F,header=F,sep='\t')
  df_normal_binned <- read.table(file.path(input_dir,paste0(sample,".normal.binned.gz")),stringsAsFactors=F,header=F,sep='\t')
  colnames(df_normal_binned) <- c("chr","start","end","ref","gc")
  colnames(df_cleaned) <- c("chr","start","end","raw","gc")
  colnames(df_partitioned_raw) <- c("chr","start","end","N","segment_id")
  df_partitioned_raw$idx <- paste0(df_partitioned_raw$chr,df_partitioned_raw$start,df_partitioned_raw$end)
  df_cleaned$idx <- paste0(df_cleaned$chr,df_cleaned$start,df_cleaned$end)
  df_normal_binned$N <- NA
  df_normal_binned$raw <- NA
  subset_idx <- paste0(df_normal_binned$chr,df_normal_binned$start,df_normal_binned$end) %in% df_partitioned_raw$idx
  df_normal_binned$N[subset_idx] <- df_partitioned_raw$N 
  subset_idx <- paste0(df_normal_binned$chr,df_normal_binned$start,df_normal_binned$end) %in% df_cleaned$idx
  df_normal_binned$raw[subset_idx] <- df_cleaned$raw
  df_normal_binned$chr <- paste0("chr",df_normal_binned$chr)
  df_normal_binned$cor <- df_normal_binned$raw
  df_normal_binned$bin_size <- df_normal_binned$end-df_normal_binned$start

  gr_cn <- makeGRangesFromDataFrame(df_normal_binned)
  
  bin_wm <- function(j,v1,v2,idx_first,idx_last) { 
  ifelse(!is.na(idx_first[j]) & (!is.na(idx_last[j])),
         weighted.mean(v1[idx_first[j]:idx_last[j]],v2[idx_first[j]:idx_last[j]],na.rm=T),NA)
  }
  bin_med <- function(j,v,idx_first,idx_last) { 
    ifelse(!is.na(idx_first[j]) & (!is.na(idx_last[j])),
           median(v[idx_first[j]:idx_last[j]],na.rm=T),NA)
  }
  bin_best <- function(j,v1,v2,idx_first,idx_last) { 
    ifelse(!is.na(idx_first[j]) & (!is.na(idx_last[j])),
           v1[idx_first[j]:idx_last[j]][which.max(v2[idx_first[j]:idx_last[j]])], NA)
  }
  bin_sum <- function(j,v,idx_first,idx_last) { 
    ifelse(!is.na(idx_first[j]) & (!is.na(idx_last[j])),
           sum(v[idx_first[j]:idx_last[j]],na.rm=T),NA)
  }
  
  bin_resolutions <- c("wg_coarse","wg_medium","wg_fine","targeted")
  for (resolution in bin_resolutions) {
    print(resolution)
    bin_file <- paste0("hg19_",resolution,".bed")
    gr_bins <- readBed(file.path(annotations_dir,bin_file))
    idx_first <- findOverlaps(gr_bins,gr_cn,select="first")
    idx_last <- findOverlaps(gr_bins,gr_cn,select="last")
    df <- as.data.frame(gr_bins)
    colnames(df)[1] <- "chr"
    df$score <- NULL
    df$raw <- sapply(1:nrow(df),bin_sum,df_normal_binned$raw,idx_first,idx_last)
    df$cor <- df$raw
    df$N <- sapply(1:nrow(df),FUN=bin_med,df_normal_binned$N,idx_first,idx_last)
    df$N <- df$N*2/median(df$N,na.rm=T)
    df$Ns <- round(df$N)
    assign(paste0("df_",resolution),df)
  }
  if(IMPORT_CNV_VCF) {
    vcf <- readVcf(file.path(input_dir,paste0(sample,".vcf.gz")))
    exp_vcf <- expand(vcf)
    df_vcf <- as.data.frame(info(exp_vcf))
    list_seg <- strsplit(rownames(df_vcf),"[:-]")
    df_seg <- as.data.frame(matrix(unlist(list_seg),ncol=5,byrow=T))
    colnames(df_seg) <- c("Caller","GainLoss","chr","start","end")
    df_seg$N <- 2
    df_seg$N[df_cnv$GainLoss=="Gain"] <- 3
    df_seg$N[df_cnv$GainLoss=="Loss"] <- 1
    df_seg$Caller <- NULL
    df_seg$GainLoss <- NULL
    n_sample_files <- 5
  } else {   n_sample_files <- 4;  }
  # Add to index df and write sample files
  df_index_entry <- data.frame(path=c("wg_coarse","wg_medium","wg_fine","targeted","seg"),	
                         file=c("counts","counts","counts","counts","deletions"),	
                         sample=sample,	
                         resolution=c("wg_coarse","wg_medium","wg_fine","targeted","wg_medium"),	
                         is_targeted=c(F,F,F,T,F),	
                         is_wgs=c(T,T,T,F,F))
  dfs <- paste0("df_",df_index_entry$path)
  df_index_entry$path <- paste0(sample,"_",df_index_entry$path,".tsv")
  for (j in 1:n_sample_files) {
    print(df_index_entry$path[j])
    write.table(get(dfs[j]),file.path(output_dir,df_index_entry$path[j]),sep='\t',quote=F,row.names = F)
  }
  df_index <- rbind(df_index,df_index_entry[1:n_sample_files,])
}
write.table(df_index,index_file,sep='\t',quote=F,row.names = F)
