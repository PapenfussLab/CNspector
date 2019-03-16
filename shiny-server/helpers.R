library(shiny)
library(futile.logger)
library(plyr)
library(data.table)
library(R.utils)
library(RColorBrewer)

# TODO: 
# ? help switch prints proper use from a file (if it exists)
# print library versions in use
# print more diagnostics reading files when verbose/debug=T
# change verbose/debug so that debug will print stuff and update per-message but verbose does not.
#
USE_ALT_READ_TABLE <- TRUE # Replace read.table() with fread() (faster but sometimes gets format wrong)
HEAD_N <- 1000 # number of lines to read in to try and guess file format

# Application-specific types that we don't need to guess
INTEGER_COLS <- c("start","end","Start","End","start1","Start1","start2","Start2",
                  "end1","End1","end2","End2","alt_reads","ref_reads","read_depth",
                  "bins","Length","binSize","segments","refpoolsize"
                  )
CHARACTER_COLS <- c("chr","Chr","chr2","Chr2","Sample","SampleName",
                    "strand","name","bam","bed","refname","sample_type","path","file",
                    "sample","panel","resolution","refname",
                    "sex","Sex","tumour_normal"  
                    )
LOGICAL_COLS <- c("IsCellLine","IsRnaSeq","LibraryFlag","targeted","blackListed",
                  "is_wgs","whiteListed","outlier","outlier_sd","outlier_z","is_targeted"
                  )

# -------------------------------------------------------------------
# Read bed-like or other tsv tables tables with read.table()-like behaviour
# Use some pre-defined column types to avoid making the wrong assumptions about type.
# -------------------------------------------------------------------
read.table.bed <- function (file,...){
  if(USE_ALT_READ_TABLE==FALSE) { return (read.table(file,...)) }
  dots <- list( file, ... )
  if( ("header" %in% names(dots)) && (dots[["header"]]==FALSE) ){ return(read.table(file,...)) }
  df_head <- read.table(file,nrows=HEAD_N,...)
  if( "colClasses" %in% names(dots) ){
    colClasses <- dots[["colClasses"]]
    # maybe override? not for now
  } else {
    colClasses <- character()
    for (colname in colnames(df_head)) {
      colClasses[[colname]] <- class(df_head[[colname]])
      if( colname %in% CHARACTER_COLS ){ colClasses[[colname]] <- "character" }
      else if( colname %in% LOGICAL_COLS ){ colClasses[[colname]] <- "logical" }
      # else if( colname %in% INTEGER_COLS ){ colClasses[[colname]] <- "integer" }
    }
    dots[["colClasses"]] <- colClasses
  }
  dots[["data.table"]] <- FALSE
  df <- do.call("fread", dots)
  return(df)
}

# -------------------------------------------------------------------
# A reporting function that can write to the screen if needed
# -------------------------------------------------------------------
cnb_echo <- function(gv,output,func,str,...){
  if(!exists("isolate")) { isolate <- function (f) { return (f)} }
  if(gv$init_tab$verbose || gv$init_tab$debug) { gv$log <- paste0(isolate(gv$log),func(str,...)) }
  if(gv$init_tab$debug && gv$haveShiny) { output$error_message_text <- renderText({ gv$log }) }
  return(gv) 
}  

# -------------------------------------------------------------------
# Initialise globals that either don't change or are defaults and don't change after init_globals()
# -------------------------------------------------------------------
init_constant_globals <- function(...){
  args <- list(...)
  gv <- list() # global variables
  gv$init_table <- list()
  gv$init_table$infile_path <- NULL
  gv$cyto.anno <- NULL
  gv$chrom.info <- NULL
  gv$log <- ""
  gv$samples  <- c("")
  gv$exons <- NULL
  gv$gene_names <- NULL
  gv$ok <- FALSE
  gv$haveShiny <- exists("renderText") # a function in shiny used for error reporting
  # Warning - FLAGS are bitwise but they also can be an index to gv$resolutions and gv$file_types
  gv$COARSE_IDX = 1
  gv$MEDIUM_IDX = 2
  gv$FINE_IDX = 3
  gv$WG_MASK = 3
  gv$COUNTS_MASK = 15
  gv$FILE_TYPE_MASK = 240
  gv$WG_FLAG = 1
  gv$TARGETED_FLAG = 4
  gv$AUX_TARGET_FLAG = 8
  gv$BAF_FLAG = 16
  gv$FUSION_FLAG = 32
  gv$SEGMENTATION_FLAG = 64
  gv$ERROR_BARS_FLAG = 128
  # Some constants
  gv$cn_different_colours <- 9
  gv$cn_max_call <- 1e4
  gv$cn_colour_scale <- list()
  gv$fusion_colour_scale <- colorRampPalette(brewer.pal(9,"Blues")[c(2,9)])(20)
  gv$fusion_colour <- brewer.pal(8,"Blues")[8]
  gv$cn_colour_scale <- c("black","red","grey50","blue","green","yellowgreen","magenta","brown",rep("orange",gv$cn_max_call-gv$cn_different_colours+2))
  names(gv$cn_colour_scale) <- as.factor(0:(length(gv$cn_colour_scale)-1))
  gv$plot2_colour <- "grey80"
  gv$plot3_colour <- "grey80"
  gv$segmentation_colour <- "black"
  gv$text_offset_fraction <- 0.1
  gv$fusion_line_colour <- "darkviolet"
  gv$fusion_line_size <- 0.5
  gv$cn_min_val <- 0
  gv$resolutions <- c("wg_coarse","wg_medium","wg_fine","targeted","","","","aux_target")
  gv$file_type[1:16] <- "counts"
  gv$file_type[gv$BAF_FLAG] <- "baf"
  gv$file_type[gv$FUSION_FLAG] <- "fusion"
  gv$file_type[gv$SEGMENTATION_FLAG] <- "deletions"
  gv$data_type_names <- c("Whole Genome","Whole Genome","Whole Genome","Targeted","","","","Off Target")
  gv$first_index = 1 # Offset of first base in chr - 1 based index
  gv$gene_names_filename <- NULL
  gv$exons_filename <- NULL
  gv$targeted_gene_names_filename <- NULL
  gv$yAutosomeRatioThreshold  <- 0.125
  gv$xAutosomeRatioThreshold <- 0.75
  gv$hidden_points <- data.frame(Ns = 0:(gv$cn_different_colours-1),N = 0,x = -1,x1 = -1,x2 = -1,Nmin=0,Nmax=4,type=0)
  gv$hidden_points <- NULL
  gv$init_table$plot1_point_size <- 1
  gv$init_table$plot2_point_size <- 0.1
  gv$init_table$plot3_point_size <- 0.1
  gv$init_table$cn_max_display_val <- 4
  gv$init_table$cn_min_display_val <- -0.2
  gv$init_table$cn_min_display_val_facted <- -1.2
  gv$init_table$have_wg <- FALSE
  gv$init_table$have_targeted <- FALSE
  gv$init_table$have_aux_target <- FALSE
  gv$init_table$have_segmentation <- FALSE
  gv$init_table$have_fusion <- FALSE
  gv$init_table$have_baf <- FALSE
  gv$init_table$have_error_bar <- FALSE
  gv$init_table$locked <- FALSE
  gv$init_table$max_fusions <- 50
  gv$init_table$max_samples <- 100
  gv$init_table$display_sample <- 1 # or can be sample name
  gv$init_table$sample <- NULL
  gv$init_table$run <- NULL  
  gv$init_table$reference_path <- ""  
  gv$init_table$base_path <- ""
  gv$init_table$collection_root_path <- "/pathology/NGS/Samples" # site specific path to all collections
  gv$init_table$collection <- NULL
  gv$init_table$panel <- NULL  
  gv$init_table$min_probe_coverage <- 1
  gv$init_table$chr <- "chr17"
  gv$init_table$chr_label <- "17" # the label that is printed on the plot
  gv$init_table$start <- 7571720
  gv$init_table$end <- 7590868
  gv$init_table$arm <- "p"  
  gv$init_table$display_gene <- "TP53"  
  gv$init_table$gc_correct_wg <- FALSE
  gv$init_table$gc_correct_targeted <- FALSE
  gv$init_table$gc_correct_aux_target <- FALSE
  gv$init_table$use_reference_wg <- FALSE
  gv$init_table$use_reference_targeted <- FALSE
  gv$init_table$use_reference_aux_target <- FALSE
  gv$init_table$display_all_gene_names = FALSE
  gv$init_table$display_wg_data = TRUE
  gv$init_table$display_cytobands = TRUE
  gv$init_table$display_aux_target_data = TRUE
  gv$init_table$display_targeted_data = TRUE
  gv$init_table$display_baf <- FALSE
  gv$init_table$display_fusion <- FALSE
  gv$init_table$display_segmentation <- FALSE
  gv$init_table$is_gene <- FALSE
  gv$init_table$rna <- FALSE
  gv$init_table$rna_coverage_cutoff <- 10
  gv$init_table$rna_normalisation_bins <- 1000
  gv$init_table$use_pre_built_annotations <- FALSE
  gv$init_table$multi_sample_mode <- FALSE
  gv$init_table$index_path <- ""
  gv$init_table$n_scale <- 1.0
  gv$init_table$debug <- FALSE
  gv$init_table$verbose <- FALSE
  #gv$init_table$infile_path <- "default.tsv"
  #gv$init_table$verbose <- T
  #gv$init_table$debug <- T
  gv[names(args)] <- args[names(args)]
  return(gv)
}

# -------------------------------------------------------------------
# Initialise globals from URL
# -------------------------------------------------------------------
init_globals <- function(session_info,session,output){
# For stand-alone debugging outside of the Shiny framework
  if(!exists("session")) {session <- NULL} 
  if(!exists("session_info")) { session_info <- NULL }
  if(!exists("output")) { output <- list(error_message="") }    
  if(!exists("isolate")) { isolate <- function (f) { return (f)} }
  if(!exists("ftry")) { ftry <- try }
  gv <- init_constant_globals()
  
  if(is.null(session_info)) {
    gv$init_table$debug <- TRUE
    gv$init_table$verbose <- TRUE
    gv$init_table$log <- ""
    return(gv)
  }
   
  # Plots are number 1, 2, 3 from the top, and the table is 4
  # reactiveValues - displayOptions2, displayOptions3, displayOptions4 are named after the plot which depends on them  
  # controlsAll - controls which can change all plots 
  # Point type is a bitwise flag - which is also used to index gv$resolutions below (oops) 
  if(!is.null(session_info)) {
    idx_not_nominated_string <- !(names(session_info) %in% c("run", "sample", "display_sample", "panel", "gene_id", "base_path", "sample_path"))
    idx_logical  <- suppressWarnings(!is.na(as.logical(session_info[names(session_info)])) & idx_not_nominated_string)
    idx_numeric <- suppressWarnings(!is.na(as.numeric(session_info[names(session_info)])) & idx_not_nominated_string)
    idx_str <- suppressWarnings(!is.na(as.character(session_info[names(session_info)]))) & (!idx_logical) & (!idx_numeric)
    gv$init_table[names(session_info)[idx_logical]] <- as.logical(session_info[names(session_info)[idx_logical]])
    gv$init_table[names(session_info)[idx_numeric]] <- as.numeric(session_info[names(session_info)[idx_numeric]])
    gv$init_table[names(session_info)[idx_str]] <- as.character(session_info[names(session_info)[idx_str]])
    if(sum(names(session_info)=="reference_path")==0 && sum(names(session_info)=="base_path")==1) {
      gv$init_table$reference_path <- gv$init_table$base_path
    }
  }
  
  gv <- cnb_echo(gv,output,flog.info,session_info$log)
  
  if(!is.null(gv$init_table$run) && is.null(gv$init_table$panel) && is.null(gv$init_table$collection)) {
    gv <- cnb_echo(gv,output,flog.info,"Have infile_path")
    gv$init_table$infile_path <- file.path(gv$init_table$base_path,paste0(gv$init_table$run,".tsv"))
  }

  if(!(is.null(gv$init_table$run) || is.null(gv$init_table$panel) || is.null(gv$init_table$collection) || is.null(gv$init_table$sample))) {
    gv <- cnb_echo(gv,output,flog.info,"Constructing infile_path from components")
    gv$init_table$reference_path <- ""
    gv$init_table$base_path <- ""
    gv$init_table$infile_path <- file.path(gv$init_table$base_path,
                                           gv$init_table$collection_root_path,
                                           gv$init_table$collection,
                                           gv$init_table$run,
                                           gv$init_table$sample,
                                           gv$init_table$index_path,
                                           paste0(gv$init_table$panel,"_",gv$init_table$run,"_",gv$init_table$sample,".tsv"))
  }
  gv <- cnb_echo(gv,output,flog.info,"base_path = %s",gv$init_table$base_path)
  gv <- cnb_echo(gv,output,flog.info,"reference_path = %s",gv$init_table$reference_path) 
  gv <- cnb_echo(gv,output,flog.info,"collection_root_path = %s",gv$init_table$collection_root_path)        
  gv <- cnb_echo(gv,output,flog.info,"collection = %s",gv$init_table$collection)        
  gv <- cnb_echo(gv,output,flog.info,"run = %s",gv$init_table$run)        
  gv <- cnb_echo(gv,output,flog.info,"sample = %s",gv$init_table$sample)        
  gv <- cnb_echo(gv,output,flog.info,"Reading: %s",gv$init_table$infile_path)    

  if(is.null(gv$init_table$infile_path)) {
    gv <- cnb_echo(gv,output,flog.warn,"Need to specify an input file.")
    return(gv)
  }
  # Read table on the run to display
  df_files <- ftry(read.table.bed(gv$init_table$infile_path,header = TRUE,stringsAsFactors = FALSE,sep="\t"))
  if(length(df_files)<=1) {
    gv <- cnb_echo(gv,output,flog.error,"Unable to read config file: %s",gv$init_table$infile_path)
    return(gv)
  } else {
    gv <- cnb_echo(gv,output,flog.info,"Read config file: %s",gv$init_table$infile_path)
  }
  
  cols <- colnames(df_files)
  df_files$sample <- as.character(df_files$sample) # some names may look like numbers - fix those
  unique_names <- paste0(df_files$file,".",df_files$sample,".",df_files$resolution)
  idx <- (df_files$file=="reference")
  unique_names[idx] <- df_files$sample[idx]
  rownames(df_files) <- unique_names

  gv <- cnb_echo(gv,output,flog.info,"colnames(df_files): %s",paste(colnames(df_files),collapse = " , "))
  gv <- cnb_echo(gv,output,flog.info,"rownames(df_files): %s",paste(rownames(df_files),collapse = " , "))
  gv <- cnb_echo(gv,output,flog.info,"df_files$path: %s",paste(df_files$path,collapse = " , "))
  gv <- cnb_echo(gv,output,flog.info,"df_files$sample: %s",paste(df_files$sample,collapse = " , "))
  
  chrom.info.file <- file.path(gv$init_table$reference_path,df_files$path[grepl("chrom",df_files$sample)])
  if(!file.exists(chrom.info.file)) {
    gv <- cnb_echo(gv,output,flog.error,"Unable to read chrom info file: %s",chrom.info.file)
    return(gv)
  }
  chrom.info <- read.table.bed(chrom.info.file,header=FALSE,stringsAsFactors = FALSE)
  colnames(chrom.info) <- c("chr","size","total")
  rownames(chrom.info) <- chrom.info$chr
  chr_idx <- c(paste0("chr",1:22),"chrX","chrY")
  chrom.info <- chrom.info[chr_idx,]
  cs <- cumsum(as.numeric(chrom.info$size))
  num_rows <- length(cs)
  chrom.info$total[2:num_rows] <- cs[1:(num_rows-1)]
  chrom.info$total[1] <- 0
  chrom.info$total <- as.numeric(chrom.info$total)
  cyto.text.file <- file.path(gv$init_table$reference_path,df_files$path[grepl("cyto",df_files$sample)])
  
  if(!file.exists(cyto.text.file)) {
    gv <- cnb_echo(gv,output,flog.error,"Unable to read cyto file: %s",cyto.text.file)
    return(gv)
  }
  cyto.text <- read.table.bed(cyto.text.file,header=FALSE,stringsAsFactors = FALSE)
  colnames(cyto.text) <- c("chr", "start", "end", "name", "gieStain")
  unique_samples <- unique(df_files$sample[df_files$file!="reference"])
  
  if(!is.character(gv$init_table$display_sample) && gv$init_table$display_sample < gv$init_table$max_samples) {
    gv$init_table$display_sample = unique_samples[gv$init_table$display_sample]
  }
  gv$init_table$display_sample <- as.character(gv$init_table$display_sample)
  gv <- cnb_echo(gv,output,flog.info,"unique_samples: %s",paste(unique_samples,collapse = " , "))
  gv <- cnb_echo(gv,output,flog.info,"gv$init_table$display_sample = %s",gv$init_table$display_sample)
  
  # Make a list of globals carrying all these tables, and put some display variables
  if(!gv$init_table$locked) {
    idx <- !(grepl("chrom",df_files$sample) | grepl("cyto",df_files$sample))  
  } else{
    idx <- (grepl("annotation",df_files$sample) | grepl(gv$init_table$display_sample,df_files$sample))
  }
  if(sum(idx)<=1) {
    gv <- cnb_echo(gv,output,flog.error,"Not enough matching files in inputfile. Sample col: %s",paste(df_files$sample,collapse = " , "))
    return(gv)
  }
  paths <- file.path(ifelse(df_files$file=="reference",gv$init_table$reference_path,gv$init_table$base_path),df_files$path)
  filelist <- paths[idx]
  for (filepath in filelist) {
    if(!file.exists(filepath)) {
      gv <- cnb_echo(gv,output,flog.error,"Required file not present: %s",filepath)
      return(gv)
    }
  }
  fapply_1 <- function(fn, ci) { 
    display_fn <- sub("_.*_","~",fn)
    if(!is.null(progress)) { progress$inc(amount = 1, message = "Reading ", detail = display_fn) }
    gv <- cnb_echo(gv,output,flog.info,"Reading: %s",display_fn)
    df <- ftry(read.table.bed(fn,header=TRUE,stringsAsFactors = FALSE,sep = "\t"))
    if(length(df)==1) {
      gv <- cnb_echo(gv,output,flog.warn,"Unable to read file: %s.",display_fn)
      return(NULL)
    }
    # Only want chr1-22,X,Y
    if(sum(colnames(df)=="chr2")>0) { 
      df <- df[!grepl("_",df$chr) & !grepl("MT",df$chr) & !grepl("MT",df$chr2)  & !grepl("GL",df$chr)  & !grepl("GL",df$chr2),]
    } else {
      df <- df[!grepl("_",df$chr) & !grepl("MT",df$chr)  & !grepl("GL",df$chr),]
    }
    
    # Add display offset into genomic co-ords
    chr_offset1 <- ci[as.character(df$chr),"total"]
    if(sum(colnames(df)=="chr2")>0) { 
      chr_offset2 <- ci[as.character(df$chr2),"total"] 
    } else { 
      chr_offset2 <- chr_offset1 
    }
    df$x1 <- df$start + chr_offset1 
    df$x2 <- df$end   + chr_offset2
    df$x <- 0.5 * (df$x1 + df$x2) 
    return(df)
  }
  progress <- NULL
  if(!is.null(session)) { progress <- shiny::Progress$new(session, min=1, max=length(filelist)) }
  if(!is.null(session)) { progress$set(message = 'Reading data', detail = 'Please wait...')  }
  m <- ftry(lapply(filelist,fapply_1,chrom.info))
  if(length(m)!=length(filelist)) {
    gv <- cnb_echo(gv,output,flog.error,"Broke while reading inputs df_files on or after %s",file_list[length(m)])
    return(gv)
  }
  names(m) <- rownames(df_files)[idx]
  gv$init_table$reference_samples <- unique_samples
  if(length(unique_samples)>1) { gv$init_table$multi_sample_mode <- T }
  
  # Make a list of globals carrying all these tables, and put in display "x" variable
  gv$samples  <- unique_samples
  gv$m <- m
  gv$filelist <- df_files
  
  # chromosome arm info
  if(is.null(cyto.text)) {
    cytoFile <- df_files$path[grepl("cyto",df_files$file)]
    cytoFile <- file.path(gv$init_table$reference_path,cytoFile)
    cyto.text <- ftry(read.table.bed(cytoFile,header=TRUE,stringsAsFactors = FALSE))
    if(length(cyto.text)==1) {
      gv <- cnb_echo(gv,output,flog.error,"Unable to read cytobands file: %s",cytoFile)
      return(gv)
    }
  }
  cyto.pq <- cyto.text[cyto.text$gieStain=="acen",]
  idx <- grepl("p",cyto.pq$name)
  cyto.pq[idx,"start"] <- gv$first_index
  idx <- grepl("q",cyto.pq$name)
  chr <- cyto.pq[idx,"chr"]
  chr_end <- chrom.info[chr,"size"] 
  cyto.pq[idx,"end"] <- chr_end
  cyto.pq$seqnames <- factor(cyto.pq$chr)
  cyto.pq$N <- 2
  cyto.pq$Ns <- 2
  cyto.anno <- cyto.pq[grepl("p",cyto.pq$name),]
  row.names(cyto.anno) <- cyto.anno$chr
  cyto.anno$N <- gv$init_table$cn_max_display_val
  cyto.anno$chr_label <- gsub("chr", "", as.character(cyto.anno$chr))
  cyto.anno$start_p <- chrom.info[as.character(cyto.anno$chr),"total"]
  cyto.anno$end_p   <- cyto.anno[cyto.anno$chr,"end"]   + chrom.info[as.character(cyto.anno$chr),"total"]
  cyto.anno$start_q <- cyto.anno[cyto.anno$chr,"end"]   + chrom.info[as.character(cyto.anno$chr),"total"]
  cyto.anno$end_q   <- chrom.info[as.character(cyto.anno$chr),"size"]  + chrom.info[as.character(cyto.anno$chr),"total"]
  cyto.anno$x <- cyto.anno$end_p
  cyto.anno <- cyto.anno[chr_idx,]
  gv$cyto.anno <- cyto.anno
  gv$chrom.info <- chrom.info
  gv <- cnb_echo(gv,output,flog.info,"Filling exon annotation table")
  if(is.null(gv$gene_names_filename)) {gv$gene_names_filename <- file.path(dirname(as.character(df_files[rownames(df_files)=="wg_annotation_file","path"])),"prebuilt_gene_names_table.tsv")}
  if(is.null(gv$exons_filename)) {gv$exons_filename <- file.path(dirname(as.character(df_files[rownames(df_files)=="wg_annotation_file","path"])),"prebuilt_exons_table.tsv")}
  if(is.null(gv$targeted_gene_names_filename)) {gv$targeted_gene_names_filename <- file.path(dirname(as.character(df_files[rownames(df_files)=="targeted_annotation_file","path"])),"prebuilt_targeted_gene_names.tsv")}
  if(gv$init_table$use_pre_built_annotations && 
      file.exists(gv$gene_names_filename) && file.exists(gv$exons_filename) && 
      file.exists(gv$targeted_gene_names_filename)) {
    gv$gene_names <- read.table.bed(gv$gene_names_filename,header = TRUE)    
    gv$exons <- read.table.bed(gv$exons_filename,header = TRUE)    
    gv$targeted_gene_names <- read.table.bed(gv$targeted_gene_names_filename,header = TRUE)    
  } else {
    if(!is.null(progress)) { progress$inc(amount = 1, message = "Preparing annotations - ", detail = "exons") }
    if(is.null(gv$m$ref.targeted_exons)) {gv$targeted_exons <-gv$m$targeted_annotation_file } else { gv$targeted_exons <- gv$m$ref.targeted_exons}
    gv$targeted_exons$targeted <- TRUE
    if(is.null(gv$m$ref.all_exons)) {gv$exons <- gv$m$wg_annotation_file} else{ gv$exons <- gv$m$ref.all_exons}
    gv$exons$targeted <- FALSE
    idx_targeted <- gv$exons$name %in% gv$targeted_exons$name
    gv$exons$targeted[idx_targeted] <- TRUE
    gv$exons$N <- 0
    gv$gene_names <- gv$exons[gv$exons$exon_number==1,]
    gv$gene_names$N <- gv$gene_names$N + 0.5
    gv$exons$len <-  gv$exons$end - gv$exons$start
    gv$targeted_gene_names <- gv$gene_names[gv$gene_names$targeted==TRUE,]
  }
  # Make cytoband df for rendering in plot 2
  gv <- cnb_echo(gv,output,flog.info,"Preparing annotations.")
  if(!is.null(progress)) { progress$inc(amount = 1, message = "Preparing annotations - ", detail = "cytobands") }
  cyto.text$len <- cyto.text$end - cyto.text$start + 1 
  cyto.text$x <- cyto.text$start +  chrom.info[as.character(cyto.text$chr),"total"] 
  cyto.text$x1 <- cyto.text$x
  cyto.text$x2 <- cyto.text$end +  chrom.info[as.character(cyto.text$chr),"total"] 
  cyto.text$colour <- "white" # gneg
  cyto.text[grepl("gpos",cyto.text$gieStain),"colour"] <- "grey50" 
  cyto.text[grepl("gpos25",cyto.text$gieStain),"colour"] <- "grey90" 
  cyto.text[grepl("gpos50",cyto.text$gieStain),"colour"] <- "grey50" 
  cyto.text[grepl("gpos75",cyto.text$gieStain),"colour"] <- "grey20" 
  cyto.text[grepl("gpos100",cyto.text$gieStain),"colour"] <- "black" 
  cyto.text[grepl("gvar",cyto.text$gieStain),"colour"] <- "indianred" 
  cyto.text[grepl("stalk",cyto.text$gieStain),"colour"] <- "skyblue1" 
  cyto.text[grepl("acen",cyto.text$gieStain),"colour"] <- "red"
  cyto.text$name <- gsub("p","", cyto.text$name)
  cyto.text$name <- gsub("q","", cyto.text$name)
  cyto.text$N <- 0
  gv$cyto.text <- cyto.text
  
  # Add display offset into genomic co-ords
  display_offset <- gv$chrom.info[as.character(gv$init_table$chr),"total"]
  gv$init_table$start <- gv$init_table$start + display_offset
  gv$init_table$end <- gv$init_table$end + display_offset
  # Far right of plot 1
  gv$init_table$max_x <- max(gv$cyto.anno[,"end_q"])

  # Once chr and exon info has been read in then chromosome arm and gene co-ords can be found
  if(!is.null(session_info$gene_id)) {
    initial_exons <- gv$exons[(gv$exons$name==session_info$gene_id),]
    if(nrow(initial_exons) > 0) {
      # remove duplicate gene names - just use the first
      gv$init_table$chr  <- as.character(initial_exons[1,"chr"])
      initial_exons <- initial_exons[(initial_exons$chr==gv$init_table$chr),]
      # co-ords come from these exons
      gv$init_table$start <- min(initial_exons$x - initial_exons$len)
      gv$init_table$end   <- max(initial_exons$x + initial_exons$len)
      gv$init_table$is_gene <- TRUE
      gv$init_table$display_gene_id <- session_info$gene_id
    }
  }
  if(gv$init_table$start > gv$cyto.anno[as.character(gv$init_table$chr),"end_p"]) { 
    gv$init_table$arm <- "q" 
    gv$init_table$start_arm <- gv$cyto.anno[as.character(gv$init_table$chr),"start_q"] 
    gv$init_table$end_arm   <- gv$cyto.anno[as.character(gv$init_table$chr),"end_q"] 
  } else {
    gv$init_table$arm <- "p" 
    gv$init_table$start_arm <- gv$cyto.anno[as.character(gv$init_table$chr),"start_p"] 
    gv$init_table$end_arm   <- gv$cyto.anno[as.character(gv$init_table$chr),"end_p"] 
  }
  gv$init_table$have_wg <- sum(grepl(gv$resolutions[gv$WG_FLAG],gv$filelist$resolution) & gv$filelist$file=="counts")>0
  gv$init_table$have_targeted <- sum(grepl(gv$resolutions[gv$TARGETED_FLAG],gv$filelist$resolution) & gv$filelist$file=="counts")>0
  gv$init_table$have_aux_target <- sum(grepl(gv$resolutions[gv$AUX_TARGET_FLAG],gv$filelist$resolution) & gv$filelist$file=="counts")>0
  gv$init_table$have_baf <- sum(gv$filelist$file=="baf")>0
  gv$init_table$have_segmentation <- sum(gv$filelist$file=="deletions")>0
  gv$init_table$have_fusion <- sum(gv$filelist$file=="fusion")>0
  if(!is.null(session)) { progress$close() }
  gv$ok <- TRUE
  return(gv)
}

