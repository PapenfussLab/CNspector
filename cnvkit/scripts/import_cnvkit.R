#!/usr/bin/Rscript --vanilla 
# -------------------------------------------------------------------
# Converts cnvkit output into tables displayable by CNspector
# -------------------------------------------------------------------
library(genomation)
library(GenomicRanges)
library(stats)

annotations_dir <- "../annotations/preprocessing"
input_dir <- "cnvkit_output"
output_dir <- "CNspector"
if(!dir.exists(output_dir)) { dir.create(output_dir,recursive=T) }

print(file.path(annotations_dir))
index_file <- file.path(output_dir,"index.tsv")
df_index <- data.frame(path=c("wg_transcripts.bed","cytoBand.txt","chromInfo.txt","targeted_transcripts.bed"),	
                       file="reference",	
                       sample=c("wg_annotation_file","cytobands_file","chrom_info_file","targeted_annotation_file"),	
                       resolution=c("wg","wg","wg","targeted"),	
                       is_targeted=c(F,F,F,T),	
                       is_wgs=c(T,T,T,F))
sample_names <- unique(gsub("\\..*","",dir(path=input_dir,pattern="cnr$")))
for (sample in sample_names) {
  print(sample)
  df_call_cns <- read.table(file.path(input_dir,paste0(sample,".call.cns")),stringsAsFactors=F,header=T,sep='\t')
  df_cnr <- read.table(file.path(input_dir,paste0(sample,".cnr")),stringsAsFactors=F,header=T,sep='\t')
  colnames(df_cnr) <- c("chr","start","end","gene","raw","log2","weight")
  df_cnr$gene <- NULL
  df_cnr$chr <- paste0("chr",df_cnr$chr)
  df_cnr$N <- 2^(1+df_cnr$log2)
  df_cnr$cor <- df_cnr$raw
  colnames(df_call_cns) <- c("chr","start","end","gene","log2","freq","N","cn1","cn2","raw","bins","weight")
  df_call_cns$chr <- paste0("chr",df_call_cns$chr)
  
  gr_cn <- makeGRangesFromDataFrame(df_cnr)
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
    df$score=NULL
    colnames(df)[1] <- "chr"
    df$raw <- sapply(1:nrow(df),bin_sum,df_cnr$raw,idx_first,idx_last)
    df$cor <- df$raw
    df$N <- sapply(1:nrow(df),FUN=bin_best,df_cnr$N,df_cnr$weight,idx_first,idx_last)
    df$N <- df$N*2/median(df$N,na.rm=T)
    df$weight <- sapply(1:nrow(df),FUN=bin_best,df_cnr$weight,df_cnr$weight,idx_first,idx_last)
    df$Ns <- round(df$N)
    assign(paste0("df_",resolution),df)
  }
  
  df_baf <- df_call_cns[,c("chr","start","end","gene","log2","freq","raw","bins","weight")]
  df_baf$read_depth <- df_baf$raw
  df_baf$ref_reads <- df_baf$freq * df_baf$raw
  df_baf$alt_reads <- (1-df_baf$freq) * df_baf$raw
  df_seg <- df_call_cns

  df_index_entry <- data.frame(path=c("wg_coarse","wg_medium","wg_fine","targeted","baf","seg"),	
                         file=c("counts","counts","counts","counts","baf","deletions"),	
                         sample=sample,	
                         resolution=c("wg_coarse","wg_medium","wg_fine","targeted","targeted","wg_medium"),	
                         is_targeted=c(F,F,F,T,T,F),	
                         is_wgs=c(T,T,T,F,T,F))
  dfs <- paste0("df_",df_index_entry$path)
  names(dfs) <- df_index_entry$path
  rownames(df_index_entry) <- df_index_entry$path
  df_index_entry$path <- paste0(sample,"_",df_index_entry$path,".tsv")
  # Add to index df and write sample file for this resolution
  for (j in 1:nrow(df_index_entry)) {
    print(df_index_entry$path[j])
    write.table(get(dfs[j]),file.path(output_dir,df_index_entry$path[j]),sep='\t',quote=F,row.names = F)
  }
  df_index <- rbind(df_index,df_index_entry)
}
write.table(df_index,index_file,sep='\t',quote=F,row.names = F)
