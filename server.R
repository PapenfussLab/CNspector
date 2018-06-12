library(shiny)
library(futile.logger)
library(plyr)
library(RColorBrewer)
source("helpers.R")

shinyServer(function(input, output, session) {
  readSessionInfo <- reactive({
    init_table <- list(log = "")
    cdata <- session$clientData
    cnames <- names(cdata)
    if(nchar(cdata$url_search)>2) {
      query <- parseQueryString(cdata$url_search)
      init_table$log <- paste0(init_table$log,"\nURL: ",cdata$url_search)
      init_table$log <- paste0(init_table$log,"\nParsed Query: ",paste(names(query), query, sep = " = ", collapse=", "))
      if(!is.null(query$gene_id)) {
        init_table$gene_id <- query$gene_id
      }
      if(!is.null(query$locus)) {
        locus <- query$locus
        parts <- unlist(strsplit(locus,":"))
        coords <- unlist(strsplit(parts[2],"-"))
        init_table$start <- coords[1]
        init_table$end <- coords[2]
        init_table$chr <- parts[1]
        init_table$chr_label <- init_table$chr
      }
      init_table$infile_path <- query$file
      init_table$display_sample <- query$sample
      init_table[names(query)] <- query[names(query)]
    }
    return(init_table)
  })
  flog.info("isolate(readSessionInfo())")
  session_info <- isolate(readSessionInfo())

  # Make a place for output text
  output$error_message <- renderUI({ verbatimTextOutput("error_message_text") })  
  output$error_message_text <- renderText({ gv$log })
  
  # Read in files and make tables
  gv <- tryCatch(
    init_globals(session_info,session,output),
    error=function(cond) { 
      message("Initialisation failed.")
      gv <- init_constant_globals()
      output$error_message_text <- renderText({ gv$log })
      flog.error("tryCatch")
      return(gv)
    }
  )
  # On failure, say what happened
  if(gv$ok==FALSE) {
    gv <- cnb_echo(gv,output,flog.error,"Broke during initialisation.")
    gv <- cnb_echo(gv,output,flog.warn,"Usage: http://shiny.server.url/cnb/?path/to/file=input_file.tsv")  
    gv <- cnb_echo(gv,output,flog.warn,"This can be followed by several 'option=val' pairs separated by & (no spaces though). For example") 
    gv <- cnb_echo(gv,output,flog.warn,"http...cnb/?file=input_file.tsv&sample=12K439&gene_id=BRCA1")
    gv <- cnb_echo(gv,output,flog.warn,"or")
    gv <- cnb_echo(gv,output,flog.warn,"http...cnb/?file=input_file.tsv&sample=12K439&locus=3:3000000-4000000&debug=TRUE")
    gv <- cnb_echo(gv,output,flog.warn,"Supported options:")
    gv <- cnb_echo(gv,output,flog.warn,"reference_path=REFERENCE_PATH")
    gv <- cnb_echo(gv,output,flog.warn,"base_path=BASE_PATH")
    gv <- cnb_echo(gv,output,flog.warn,"run=RUN_NAME")
    gv <- cnb_echo(gv,output,flog.warn,"sample=SAMPLE_NAME")
    gv <- cnb_echo(gv,output,flog.warn,"gene_id=GENE_NAME")
    gv <- cnb_echo(gv,output,flog.warn,"locus=CHR:START-END")
    gv <- cnb_echo(gv,output,flog.warn,"debug=TRUE|FALSE")
    gv <- cnb_echo(gv,output,flog.warn,"verbose=TRUE|FALSE")
    gv <- cnb_echo(gv,output,flog.warn,"locked=TRUE|FALSE")
    output$error_message_text <- renderText({ gv$log })
    return;
  } else if(gv$init_tab$debug) {
    output$error_message_text <- renderText({ gv$log })
  }
  
  progress <- shiny::Progress$new(session, min=0, max=1)
  progress$set(message = 'Preparing for display', detail = '...')
  
  # -------------------------------------------------------------------
  # Make DF with/without reference and with/without GC correction or CN calling
  # -------------------------------------------------------------------
  setalldf <- function(gv,display_sample,reference_samples,use_reference,use_gc_correction,data_type){
    # reference_samples = list of samples to build reference out of                       
    # Incoming cor, ann on freec.  Various for me.
    # Outgoing raw, cor, y/N, Ns, type, ref 
    progressdf <- shiny::Progress$new(session, min=0, max=1)
    progressdf$set(message = 'Preparing data for display', detail = '...')
    df_files <- gv$filelist 
    idx <- (gv$resolutions[data_type] == df_files$resolution | !grepl("counts",df_files$file)) & 
          display_sample == df_files$sample & 
          gv$file_type[data_type] == df_files$file
    num_matches <- sum(idx)
    file_type <- gv$file_type[data_type]
    if(use_gc_correction) { fieldName <- "cor" } else { fieldName  <- "raw" } 
    if(bitwAnd(data_type,gv$WG_MASK)) { typeFlag <- gv$WG_FLAG } else { typeFlag  <- data_type}
    
    cnb_echo(gv,output,flog.info,(paste0("data_type =  ",data_type,
                                         " typeFlag =  ",typeFlag,
                                         " res = ", gv$resolutions[bitwAnd(data_type,gv$WG_MASK)],
                                         " sample =  ",display_sample,
                                         " Found ",num_matches,
                                         " matches for ",display_sample,
                                         " and ",file_type)))
    
    if(num_matches==1) {
      plot_name <- rownames(df_files)[idx]
      df <- gv$m[[plot_name]]
      cols <- colnames(df)
      if(use_reference==FALSE) {
        # We know it's freec
        if(sum(grepl("ann",cols))==1 && sum(grepl("cor",cols))==1 && bitwAnd(data_type,gv$WG_MASK)>0) {
          colnames(df)[colnames(df)=="cor"] <- "N"
          colnames(df)[colnames(df)=="ann"] <- "Ns"
          df$type <- typeFlag
          df$raw <- 0
          df$ref <- 0
          df$cor <- 0
          progressdf$close()
          return(df)
        } else {
          if(file_type=="baf") {
            df$N <- df$freq - 1
            df$Nmin <- df$N
            df$Nmax <- df$N
            df$Ns <- round(1 + 2 * df$freq)
            
            colnames(df)[colnames(df)=="freq"] <- "BAF"
            df$raw <- df$ref_reads + df$alt_reads
            df$sd_estimate  <- sqrt(df$raw)
            df$Nsd <- df$sd_estimate
          }
          else if(file_type=="fusion") {
            df$nSupport <- as.integer(df$nSupport)
            idx <- order(df$nSupport,decreasing=TRUE)
            df <- df[idx,]
            df$blackListed <- TRUE
            df$blackListed[1:min(gv$init_table$max_fusions,nrow(df))] <- FALSE
            df$raw <- df$nSupport
            rawN <- log10(df$nSupport)
            df$N <- (rawN-min(rawN)) * 4 / (max(rawN)-min(rawN))
            df$Ns <- round(log10(df$nSupport))
          }
          else if(file_type=="deletions") {
            df$blackListed <- FALSE
            if(sum(colnames(df)=="bstat")==1) {
              bstat_on_bins <- mean(df$bstat/df$bins,na.rm = TRUE) # quick and dirty approx for display
              if(is.na(bstat_on_bins)) { bstat_on_bins=1 }
              df$raw <- ifelse(is.na(df$bstat),bstat_on_bins * df$bins,df$bstat)
            } else {
              df$raw <- df$bins
            }
            df$Ns <- round(df$N)
            df$sd_estimate  <- sqrt(df$weighted_var)
            df$Nsd <- df$sd_estimate
          } 
          else if(file_type=="counts") { 
            if(gv$init_table$rna==TRUE) {
              # total_counts_ratio <- sum(df$raw,rm.na=TRUE)/sum(df$ref_median,na.rm=TRUE)
              low_counts_idx <- df$ref_median < gv$init_table$rna_coverage_cutoff | df$raw < gv$init_table$rna_coverage_cutoff
              r <- log10(df$raw/df$ref_median)
              df$raw>gv$init_table$rna_coverage_cutoff
              h <- hist(r[!low_counts_idx],breaks=gv$init_table$rna_normalisation_bins,plot=FALSE)
              # rmax is index to the mode
              rmax <- which.max(h$counts)
              cs <- cumsum(h$density)
              # rzero is the index to 1 percentile (hopefully past noise)
              rzero <- which(cs>0.01)[1]
              df$N <- (r - h$mids[rzero]) * 2 / (h$mids[rmax]-h$mids[rzero])
              df$N[is.na(df$N)] <- 0
              # get rid of low counts bins
              df$N[low_counts_idx] <- 0
		#   gv <- cnb_echo(gv,output,flog.warn,"locked=TRUE|FALSE")
		#    output$error_message_text <- renderText({ gv$log })

              gv <- cnb_echo(gv,output,flog.info," mode = ")
              output$error_message_text <- renderText({ gv$log })
              #gv <- cnb_echo(gv,output,flog.info,"tcr = %f, mode = %d zero = %d low_counts = %d",
              # total_counts_ratio,h$mids[rmax],h$mids[rzero],sum(low_counts_idx,na.rom=TRUE))
            }
            if(sum(colnames(df)=="raw")==0) { df$raw  <- 1 }  
            if(sum(colnames(df)=="ref_mad")>0) { df$sd_estimate  <- (df$ref_mad/df$ref_median) * df$N }  
            if(sum(colnames(df)=="Nsd")>0) { df$sd_estimate  <- df$Nsd }  
            df$Nmax <- df$N + df$sd_estimate
            df$Nmin <- df$N - df$sd_estimate
            df$Nmin[df$Nmin<0] <- 0 
            df$Ns <- round(df$N)
          }

          if(sum(colnames(df)=="ref_mad")>0) { df$Nsd <- df$N * (df$ref_mad/df$ref_median) } # mad = 1.4826*sd so this is about sd
          if(sum(colnames(df)=="Nsd")>0) { df$Nlow <- ifelse(df$N>df$Nsd,df$N-df$Nsd,0)  }
          if(sum(colnames(df)=="Nsd")>0) { df$Nhigh<- df$N+df$Nsd }
          df$type <- typeFlag
          idx_skip <- is.na(df$N) | is.na(df$Ns)
          if(sum(colnames(df)=="blackListed")>0) { idx_skip <- idx_skip | df$blackListed }
          df <- df[!idx_skip,]
          progressdf$close()
          return(df)
        }
      }

      num_reference_samples <- length(reference_samples)
      if(use_reference==TRUE && num_reference_samples>0 &&  bitwAnd(data_type,gv$COUNTS_MASK)>0) {
        idx <- gv$resolutions[data_type] == df_files$resolution & df_files$file=="counts" & df_files$sample %in% reference_samples
        if(sum(idx)>0) {
          # refnames <- samples[idx]
          refnames <- rownames(df_files)[idx]
          is_x <- gv$m[[refnames[1]]]$chr=="chrX"
          is_y <- gv$m[[refnames[1]]]$chr=="chrY"
          is_autosome <- !(is_x | is_y)
          # Make reference
          fapply_2 <- function(fn) { gv$m[[fn]][,fieldName] }
          reftable <- sapply(refnames,fapply_2)
          # reftable[reftable==0]=1 # avoid divide by zero
          # median normalise coverage
          # fapply_mean_per_col <- function(n,reftable) { reftable[,n]/mean(reftable[,n]) }
          fapply_median_per_col <- function(n,m) { 
            v <- m[,n]; 
            is_non_zero <- (!is.na(v) & v>0)
            medianXCoverage <- median(v[is_non_zero & is_x])
            medianYCoverage <- median(v[is_non_zero & is_y])
            medianAutosomeCoverage <- median(v[is_non_zero & is_autosome])
            if(is.na(medianAutosomeCoverage) || is.null(medianAutosomeCoverage)) {
              v[] <- 0
              return(v)
            } else {
              if(is.na(medianYCoverage) || medianYCoverage < medianAutosomeCoverage * gv$yAutosomeRatioThreshold) {
                v[is_y] <- 0 # means it won't get used
              } else if (!is.na(medianYCoverage)){
                v[is_y] <- v[is_y] * 2 # everything is scaled by two later
              }
              if(is.na(medianXCoverage)) {
                cnb_echo(gv,output,flog.info,"%s: X=0",refnames[n])
                v[is_x] <- 0 # means it won't get used
              }
              else if(medianXCoverage < medianAutosomeCoverage * gv$xAutosomeRatioThreshold) {
                v[is_x] <- v[is_x] * 2 # everything is scaled by two later
                cnb_echo(gv,output,flog.info,"%s: X*=2",refnames[n])
              }
              return(v/medianAutosomeCoverage)
            }
          }
          progressdf$inc(amount = 0.1, message = "Pooled reference - ", detail = "normalising")
          m <- sapply(1:ncol(reftable),fapply_median_per_col,reftable)
          progressdf$inc(amount = 0.1, message = "Pooled reference - ", detail = "calculating")
          fapply_median_per_row <- function(n,m) { v <- m[n,]; median(v[!is.na(v) & v>0]) }
          # fapply_trimmed_mean_per_row <- function(n,m) {mean(m[n,],trim=0.1)}
          # At this point we have assembled all the references
          # Look at variance here:
          denominator <- sapply(1:nrow(reftable),fapply_median_per_row,m)
          progressdf$inc(amount = 0.1, message = "Pooled reference - ", detail = "aaplying")
          v <- df[,fieldName]
          numerator <- df[,fieldName]/median(v[!is.na(v) & v>0]) # median normalise sample
          ratio <- (numerator/denominator)
          df$N <- (2 * numerator/denominator)
          df$Ns <- round(df$N)
          df$type <- typeFlag
          df$ref <- denominator
          df$sd_estimate  <- (df$ref_mad/df$ref_median) * df$N
          df$Nmax <- df$N + df$sd_estimate
          df$Nmin <- df$N - df$sd_estimate
          df <- df[denominator>0 & !is.na(numerator) & !is.na(denominator),]
          progressdf$close()
          return(df)
        }
      }
    } else {
      progressdf$close()
      return (NULL)
    }
  }
  
  # -------------------------------------------------------------------
  #  Filter annotations in a selection
  # -------------------------------------------------------------------
  selectGenesInBox <- function(h, box,display_all) { 
    if(!is.null(h) && !is.null(box)) {
      idx <- h$N > box$ymin &  h$N < box$ymax & h$x > box$xmin & h$x < box$xmax
      if(display_all) {
        return(h[idx,])
      } else {
        return(h[idx & h$targeted,])
      }
    }
    return(NULL)
  }
  # -------------------------------------------------------------------
  #  Filter annotations in a selection
  # -------------------------------------------------------------------
  selectGenesInGenomicRange <- function(h, box,display_all) { 
    if(!is.null(h) && !is.null(box)) {
      idx <- h$x > box$xmin & h$x < box$xmax
      if(display_all) {
        return(h[idx,])
      } else {
        return(h[idx & h$targeted,])
      }
    }
    return(NULL)
  }  
  # -------------------------------------------------------------------
  #  Filter data in selection - possibly also on coverage
  # -------------------------------------------------------------------
  selectPointsInBox <- function(h, box,pointType,min_targeted_cov,hidden_points) { 
    if(is.null(h) || is.null(h$type)) {
      return(hidden_points)
    }
    idx <- (h$N > box$ymin) & (h$N < box$ymax) & (bitwAnd(h$type,pointType)!=0) & (h$raw >= min_targeted_cov) & overlapIdx(h,box)
    ret <- h[idx,]
    if(is.null(hidden_points)){
      return(ret)
    } else {
      new_columns <- colnames(h)[!colnames(h) %in% colnames(hidden_points)]
      hidden_points[,new_columns] <- h[1,new_columns]
      return(rbind.fill(ret,hidden_points))
    }
  }

  # -------------------------------------------------------------------
  #  Filter data in selection - possibly also on coverage, but allow overflow in y
  # -------------------------------------------------------------------
  selectPointsInGenomicRange <- function(h, box,pointType,min_targeted_cov,hidden_points) { 
    if(is.null(h) || is.null(h$type)) {
      return(hidden_points)
    }
    idx <- (bitwAnd(h$type,pointType)!=0) & (h$raw >= min_targeted_cov) & overlapIdx(h,box)
    ret <- h[idx,]
    if(is.null(hidden_points)){
      return(ret)
    }
    # if there are any extra columns in h then put them in hidden_points so they match
    new_columns <- colnames(h)[!colnames(h) %in% colnames(hidden_points)]
    hidden_points[,new_columns] <- h[1,new_columns]
    return(rbind.fill(ret,hidden_points))
  }

  # -------------------------------------------------------------------
  #  Find when ranges overlap in a way that means display is required.
  #  Possibly this should use IRanges
  # -------------------------------------------------------------------
  overlapIdx <- function(h,box) {
    idx <- 
      ((h$x1 <= box$xmin) & (h$x2 >= box$xmin)) | 
      ((h$x1 <= box$xmax) & (h$x2 >= box$xmax)) |
      ((h$x1 <= box$xmax) & (h$x2 >= box$xmin)) |
      ((h$x1 >= box$xmin) & (h$x2 <= box$xmax))
  }
    
  # -------------------------------------------------------------------
  # reactiveValues
  # -------------------------------------------------------------------
  cnb_echo(gv,output,flog.info,"Initialise Reactive Values")
  init_table <- reactiveValues(basename = gv$init_table$basename, 
                               dirname = gv$init_table$dirname,
                               display_sample = gv$init_table$display_sample,
                               chr = gv$init_table$chr,
                               start = gv$init_table$start,
                               end = gv$init_table$end,
                               arm = gv$init_table$arm,
                               first_time_flag = TRUE
  )

  controlsAll <- reactiveValues(cn_max_val = gv$init_table$cn_max_display_val, 
                                refresh_table = 0,
                                pending_updates = 0,
                                min_probe_coverage = gv$init_table$min_probe_coverage,
                                display_gene_id = gv$init_table$display_gene_id,
                                display_sample = gv$init_table$display_sample,
                                reference_samples = gv$init_table$reference_samples,
                                gc_correct_wg = gv$init_table$gc_correct_wg,
                                gc_correct_targeted = gv$init_table$gc_correct_targeted,
                                gc_correct_off_target = gv$init_table$gc_correct_off_target,
                                use_reference_wg = gv$init_table$use_reference_wg,
                                use_reference_targeted = gv$init_table$use_reference_targeted,
                                use_reference_off_target = gv$init_table$use_reference_off_target,
                                display_all_gene_names = gv$init_table$display_all_gene_names,
                                display_wg_data = gv$init_table$display_wg_data, 
                                display_cytobands = gv$init_table$display_cytobands, 
                                display_off_target_data = gv$init_table$display_off_target_data,
                                display_targeted_data = gv$init_table$display_targeted_data,
                                display_baf = gv$init_table$display_baf,
                                display_fusion =  gv$init_table$display_fusion,
                                display_segmentation = gv$init_table$display_segmentation,
                                display_data_bitflag = gv$WG_FLAG + gv$TARGETED_FLAG + gv$OFF_TARGET_FLAG
  )
  displayOptions1 <- reactiveValues(
                                    xmin = 0, 
                                    xmax = gv$init_table$max_x, 
                                    ymin = gv$init_table$cn_min_display_val, 
                                    ymax = gv$init_table$cn_max_display_val, 
                                    df = NULL
  )
  displayOptions2 <- reactiveValues(
                                    xmin = gv$init_table$start_arm, 
                                    xmax = gv$init_table$end_arm, 
                                    ymin = gv$init_table$cn_min_display_val, 
                                    ymax = gv$init_table$cn_max_display_val, 
                                    chr = as.character(gv$init_table$chr), 
                                    chr_label = as.character(gv$init_table$chr_label), 
                                    arm = gv$init_table$arm,
                                    df = NULL
  )
  displayOptions3 <- reactiveValues(brush = list(xmin = gv$init_table$start_arm,
                                                 xmax = gv$init_table$end_arm, 
                                                 ymin = gv$init_table$cn_min_display_val, 
                                                 ymax = gv$init_table$cn_max_display_val
                                                 ),
                                    xmin = gv$init_table$start, 
                                    xmax = gv$init_table$end, 
                                    ymin = gv$init_table$cn_min_display_val, 
                                    ymax = gv$init_table$cn_max_display_val, 
                                    chr = as.character(gv$init_table$chr), 
                                    chr_label = as.character(gv$init_table$chr_label), 
                                    arm = "p",
                                    is_gene = gv$init_table$is_gene,
                                    display_all_gene_names = TRUE, 
                                    gene_names = NULL,
                                    exons = NULL,
                                    df = NULL
  )
  displayOptions4 <- reactiveValues(
    df = NULL
  )
  df <- reactiveValues(
    targeteddf = NULL,
    off_targetdf = NULL,
    df1 = NULL,
    df2 = NULL,
    df3 = NULL,
    df_baf = NULL,
    df_segmentation = NULL,
    df_fusion = NULL,
    dfwg1 = NULL,
    dfwg2 = NULL,
    dfwg3 = NULL
  )
  # -------------------------------------------------------------------
  cnb_echo(gv,output,flog.info,"Using Reactive Construct to initialise dfs")
  set_all_dfs <- reactive({
    if(init_table$first_time_flag==FALSE) {
      controlsAll$display_data_bitflag <- 
        gv$WG_FLAG * as.integer(controlsAll$display_wg_data) + 
        gv$TARGETED_FLAG * as.integer(controlsAll$display_targeted_data) + 
        gv$OFF_TARGET_FLAG * as.integer(controlsAll$display_off_target_data) +
           as.integer(controlsAll$display_baf) * gv$BAF_FLAG + 
           as.integer(controlsAll$display_fusion) * gv$FUSION_FLAG + 
           as.integer(controlsAll$display_segmentation) * gv$SEGMENTATION_FLAG + 
           as.integer(controlsAll$display_error_bars) * gv$ERROR_BARS_FLAG
      
      #cnb_echo(gv,output,flog.info,paste0("set_all_dfs() : display_data_bitflag = ",controlsAll$display_data_bitflag,
      #             " refresh_table = ",controlsAll$refresh_table,
      #             " pending_updates = ",controlsAll$pending_updates,
      #             " first_time_flag = ",init_table$first_time_flag))
      
      what_to_update <- bitwAnd(controlsAll$refresh_table,controlsAll$display_data_bitflag)
      new_pending <- bitwAnd(controlsAll$refresh_table,bitwNot(controlsAll$display_data_bitflag))
      controlsAll$pending_updates <- bitwAnd(controlsAll$pending_updates,bitwNot(what_to_update))
      controlsAll$pending_updates <- bitwOr(controlsAll$pending_updates,new_pending)
      
      if(bitwAnd(what_to_update,gv$WG_FLAG)!=0) {
        df$dfwg1 <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_wg,controlsAll$gc_correct_wg,gv$COARSE_IDX)      
        df$dfwg2 <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_wg,controlsAll$gc_correct_wg,gv$MEDIUM_IDX)      
        df$dfwg3 <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_wg,controlsAll$gc_correct_wg,gv$FINE_IDX)
      }
      # TODO: Clean this up
      if(bitwAnd(what_to_update,gv$TARGETED_FLAG)!=0 ) {
        df$targeteddf <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_targeted,controlsAll$gc_correct_targeted,gv$TARGETED_FLAG)
      }
      if(bitwAnd(what_to_update,gv$OFF_TARGET_FLAG)!=0 ) {
        df$off_targetdf <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_off_target,controlsAll$gc_correct_off_target,gv$OFF_TARGET_FLAG)      
      }
      if(bitwAnd(what_to_update,gv$BAF_FLAG)!=0 ) {
        df$df_baf <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_targeted,controlsAll$gc_correct_targeted,gv$BAF_FLAG)
      }
      if(bitwAnd(what_to_update,gv$FUSION_FLAG)!=0 ) {
        df$df_fusion <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_targeted,controlsAll$gc_correct_targeted,gv$FUSION_FLAG)
      }
      if(bitwAnd(what_to_update,gv$SEGMENTATION_FLAG)!=0 ) {
        df$df_segmentation <- setalldf(gv,controlsAll$display_sample,controlsAll$reference_samples,controlsAll$use_reference_targeted,controlsAll$gc_correct_targeted,gv$SEGMENTATION_FLAG)
      }
      if(!is.null(df$targeteddf)) {
        si <- summary(df$targeteddf$raw, digits=3)
        si[7] <- sd(df$targeteddf$raw)
        names(si)[7] <- "SD raw"
        si[8] <- si[7]/si[3]
        names(si)[8] <- "CV raw" 
        si[8] <- NA
        if(median(df$targeteddf$cor,na.rm=TRUE)>0) {si[9] <- sd(df$targeteddf$cor,na.rm=TRUE)/median(df$targeteddf$cor,na.rm=TRUE)}
        names(si)[9] <- "CV corr" 
        si[10] <- NA
        if(median(df$targeteddf$N,na.rm=TRUE)>0) { si[10] <- sd(df$targeteddf$N,na.rm=TRUE)/median(df$targeteddf$N,na.rm=TRUE)}
        names(si)[10] <- "CV N" 
        si <- round(si*100)/100
        output$sample_summary <- renderPrint(si)
      }
      df_targeted_and_other <- rbind.fill(df$targeteddf,df$off_targetdf,df$df_baf,df$df_fusion,df$df_segmentation)
      df$df1 <- rbind.fill(df$dfwg1,df_targeted_and_other)
      df$df2 <- rbind.fill(df$dfwg2,df_targeted_and_other)
      df$df3 <- rbind.fill(df$dfwg3,df_targeted_and_other)
      displayOptions1$df <- selectPointsInGenomicRange(df$df1,displayOptions1,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
      displayOptions2$df <- selectPointsInGenomicRange(df$df2,displayOptions2,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
      displayOptions3$df <- selectPointsInBox(df$df3,displayOptions3,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
      displayOptions3$gene_names = selectGenesInGenomicRange(gv$gene_names,displayOptions3,controlsAll$display_all_gene_names)
      displayOptions3$exons = selectGenesInGenomicRange(gv$exons,displayOptions3,controlsAll$display_all_gene_names)
      displayOptions4$df <- NULL
    }
    controlsAll$refresh_table <- 0
  })

  # -------------------------------------------------------------------
  cnb_echo(gv,output,flog.info,"Render UI")
  # -------------------------------------------------------------------
  # Render UI
  # -------------------------------------------------------------------
  if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
    output$ui_reference_samples <- renderUI({  
      checkboxGroupInput( "reference_samples", 
                          label = "Targeted Pooled Reference Samples",
                          choices =  gv$init_table$reference_samples,
                          selected = gv$init_table$reference_samples,
                          inline = TRUE)
    })
  }
  if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
    output$ui_sample <- renderUI({  
      radioButtons("display_sample", 
                   label = "Display Sample",
                   choices = gv$samples,
                   selected = gv$init_table$display_sample,
                   inline = TRUE)
    })
  } else {
    output$ui_sample <- renderUI({
      radioButtons("display_sample",
                   label = "Display Sample",
                   choices = gv$init_table$display_sample,
                   selected = gv$init_table$display_sample,
                   inline = TRUE)
    })
  }
  output$min_probe_coverage <- renderUI({  
    sliderInput('min_probe_coverage', 
                'Minimum Displayable Read Support', 
                min=0, 
                max=500,
                value=gv$init_table$min_probe_coverage, 
                step=1, 
                round=0)
  })
  

  if(gv$init_table$have_wg) {
    output$display_wg_data <- renderUI({ checkboxInput(inputId = "display_wg_data", label = "WG",value = gv$init_table$display_wg_data) })
    if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
      output$use_reference_wg <- renderUI({ checkboxInput(inputId = "use_reference_wg", label = "Dynamic reference for WG",value = gv$init_table$use_reference_wg) })
      # output$gc_correct_wg <- renderUI({ checkboxInput(inputId = "gc_correct_wg", label = "GC correct WG",value = gv$init_table$gc_correct_wg) })
    }
  }
  if(gv$init_table$have_targeted) {
    output$display_targeted_data <- renderUI({ checkboxInput(inputId = "display_targeted_data", label = "Targeted",value = gv$init_table$display_targeted_data) })
    if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
      output$use_reference_targeted <- renderUI({ checkboxInput(inputId = "use_reference_targeted", label = "Dynamic reference for Targeted",value = gv$init_table$use_reference_targeted) })
      # output$gc_correct_targeted <- renderUI({ checkboxInput(inputId = "gc_correct_targeted", label = "GC correct Targeted",value = gv$init_table$gc_correct_targeted) })
    }
  }
  if(gv$init_table$have_off_target) {
    output$display_off_target_data <- renderUI({ checkboxInput(inputId = "display_off_target_data", label = "Aux Target",value = gv$init_table$display_off_target_data) })
    if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
      output$use_reference_off_target <- renderUI({ checkboxInput(inputId = "use_reference_off_target", label = "Dynamic reference for Aux Target",value = gv$init_table$use_reference_off_target) })
      # output$gc_correct_off_target <- renderUI({ checkboxInput(inputId = "gc_correct_off_target", label = "GC correct Off Target",value = gv$init_table$gc_correct_off_target) })
    }
  }
  if(gv$init_table$locked==FALSE && gv$init_table$multi_sample_mode==TRUE) {
    output$recalculate_reference <- renderUI({ actionButton(inputId = "recalculate_reference", label = "Recalculate Reference", icon = icon("cogs")) })
  }
  output$display_error_bars <- renderUI({ checkboxInput(inputId = "display_error_bars", label = "Display Errors",value = gv$init_table$display_error_bars) })
  if(gv$init_table$have_baf) {
    output$display_baf <- renderUI({ checkboxInput(inputId = "display_baf", label = "Display BAF",value = gv$init_table$display_baf) })
  }
  if(gv$init_table$have_fusion) {
  output$display_fusion <- renderUI({ checkboxInput(inputId = "display_fusion", label = "Display Fusions",value = gv$init_table$display_fusion) })
  }
  if(gv$init_table$have_segmentation) {
  output$display_segmentation <- renderUI({ checkboxInput(inputId = "display_segmentation", label = "Display Segmentation",value = gv$init_table$display_segmentation) })
  }
  # -------------------------------------------------------------------
  # Observers
  # -------------------------------------------------------------------
  cnb_echo(gv,output,flog.info,"ObserveEvents")
  observeEvent(input$cn_max_val, {
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$cn_max_val): = ",input$cn_max_val))
    controlsAll$cn_max_val <- input$cn_max_val
    displayOptions1$ymax <- controlsAll$cn_max_val
    displayOptions2$ymax <- controlsAll$cn_max_val
    displayOptions3$ymax <- controlsAll$cn_max_val
    # It may be that some previously windowed-out points need to be re-included
    displayOptions3$df <- selectPointsInBox(df$df3,displayOptions3,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
  })
  observeEvent(input$recalculate_reference,{
    cnb_echo(gv,output,flog.info,paste0("input: refresh_table = ",input$refresh_table))
    controlsAll$refresh_table <- 
      as.integer(controlsAll$use_reference_targeted) * gv$TARGETED_FLAG + 
      as.integer(controlsAll$use_reference_wg) * gv$WG_FLAG + 
      as.integer(controlsAll$use_reference_off_target) * gv$OFF_TARGET_FLAG
      set_all_dfs()
  })
  observeEvent(input$gene_id_find,{
    controlsAll$display_gene_id <- input$gene_id
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$gene_id_find): Looking for gene_id = ",input$gene_id))
    exons <- gv$exons[(gv$exons$name==controlsAll$display_gene_id),]
    if(nrow(exons)>0) {
      chr <- exons[1,"chr"]
      displayOptions3$exons <- exons[exons$chr==chr,]
      displayOptions2$chr <- chr
      displayOptions2$chr_label <- gsub("chr", "", chr) # regardless, the label is just the number
      # cnb_echo(gv,output,flog.info,head(displayOptions3$exons))
      displayOptions3$xmin <- min(displayOptions3$exons$x - displayOptions3$exons$len)
      displayOptions3$xmax <- max(displayOptions3$exons$x + displayOptions3$exons$len)
      displayOptions3$chr <- chr
      displayOptions3$is_gene  <- TRUE
      displayOptions3$chr_label <- displayOptions2$chr_label
      chr_arms <- gv$cyto.anno[gv$cyto.anno$chr==chr,] 
      if(displayOptions3$xmin > chr_arms$end_p) {
        displayOptions2$arm <- "q"
        displayOptions3$arm <- "q"
        displayOptions2$xmin <- chr_arms$start_q
        displayOptions2$xmax <- chr_arms$end_q
      } else {
        displayOptions2$arm <- "p"
        displayOptions3$arm <- "p" 
        displayOptions2$xmin <- chr_arms$start_p
        displayOptions2$xmax <- chr_arms$end_p
      }
      set_all_dfs()
    }
  })
  
  observeEvent(input$display_wg_data, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$display_wg_data): = ",input$display_wg_data))
    controlsAll$display_wg_data <- input$display_wg_data 
    if(bitwAnd(controlsAll$pending_updates,gv$WG_FLAG)!=0) {
      controlsAll$refresh_table <- as.integer(controlsAll$display_wg_data) * gv$WG_FLAG   
    }
    set_all_dfs()
    })
  observeEvent(input$display_targeted_data, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$display_targeted_data): = ",input$display_targeted_data))
    controlsAll$display_targeted_data <- input$display_targeted_data 
    if(bitwAnd(controlsAll$pending_updates,gv$TARGETED_FLAG)!=0) {
      controlsAll$refresh_table <- as.integer(controlsAll$display_targeted_data) * gv$TARGETED_FLAG
    }
    set_all_dfs()
    })
  observeEvent(input$display_off_target_data, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$display_off_target_data): = ",input$display_off_target_data))
    controlsAll$display_off_target_data <- input$display_off_target_data 
    if(bitwAnd(controlsAll$pending_updates,gv$OFF_TARGET_FLAG)!=0) {
      controlsAll$refresh_table <- as.integer(controlsAll$display_off_target_data) * gv$OFF_TARGET_FLAG
    }
    set_all_dfs()
    })
  observeEvent(input$display_all_gene_names, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$display_all_gene_names): = ",input$display_all_gene_names))
    controlsAll$display_all_gene_names <- input$display_all_gene_names 
    displayOptions3$gene_names = selectGenesInGenomicRange(gv$gene_names,displayOptions3,controlsAll$display_all_gene_names)
    # set_all_dfs()
  })
  observeEvent(input$display_cytobands, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$display_cytobands): = ",input$display_cytobands))
    controlsAll$display_cytobands <- input$display_cytobands 
  })
  observeEvent(input$gc_correct_wg, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$gc_correct_wg): = ",input$gc_correct_wg))
    controlsAll$gc_correct_wg <- input$gc_correct_wg 
    controlsAll$refresh_table <- gv$WG_FLAG
    set_all_dfs()
  })
  observeEvent(input$gc_correct_targeted, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$gc_correct_targeted): = ",input$gc_correct_targeted))
    controlsAll$gc_correct_targeted <- input$gc_correct_targeted 
    controlsAll$refresh_table <- gv$TARGETED_FLAG
    set_all_dfs()
  })
  observeEvent(input$gc_correct_off_target, { 
    controlsAll$gc_correct_off_target <- input$gc_correct_off_target
    controlsAll$refresh_table <- gv$OFF_TARGET_FLAG
    set_all_dfs()
  })
  observeEvent(input$use_reference_wg, { 
    controlsAll$use_reference_wg <- input$use_reference_wg 
    controlsAll$refresh_table <- gv$WG_FLAG
    set_all_dfs()
  })
  observeEvent(input$use_reference_targeted, { 
    controlsAll$use_reference_targeted <- input$use_reference_targeted 
    controlsAll$refresh_table <- gv$TARGETED_FLAG
    set_all_dfs()
  })
  observeEvent(input$use_reference_off_target, { 
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$use_reference_off_target): = ",input$use_reference_off_target))
    controlsAll$use_reference_off_target <- input$use_reference_off_target 
    controlsAll$refresh_table <- gv$OFF_TARGET_FLAG
    set_all_dfs()
  })

  observeEvent(input$min_probe_coverage, {
    min_probe_coverage <- input$min_probe_coverage
    cnb_echo(gv,output,flog.info,"input: observeEvent(input$min_probe_coverage): %f",controlsAll$min_probe_coverage)
    if(min_probe_coverage != controlsAll$min_probe_coverage){
      controlsAll$min_probe_coverage <- input$min_probe_coverage
      controlsAll$refresh_table <- gv$WG_FLAG + gv$TARGETED_FLAG + gv$OFF_TARGET_FLAG
      set_all_dfs()
    }
  })
  observeEvent(input$reference_samples, {
    cnb_echo(gv,output,flog.info,"input: observeEvent(input$reference_samples): ")
    controlsAll$reference_samples <- input$reference_samples
    cnb_echo(gv,output,flog.info,paste0("Selected Samples:",controlsAll$reference_samples))
    set_all_dfs() 
  })
  observeEvent(input$display_baf, {
    controlsAll$display_baf <- input$display_baf
    controlsAll$refresh_table <- gv$BAF_FLAG
    # NB: I am setting all ymins from the first plot with my "manual faceting"
    displayOptions1$ymin <- ifelse(controlsAll$display_baf,gv$init_table$cn_min_display_val_facted,gv$init_table$cn_min_display_val)
    displayOptions2$ymin <- displayOptions1$ymin
    displayOptions3$ymin <- displayOptions1$ymin
    set_all_dfs() 
  })
  observeEvent(input$display_error_bars, {
    controlsAll$display_error_bars <- input$display_error_bars
    controlsAll$refresh_table <- gv$ERROR_BARS_FLAG
    set_all_dfs() 
  })
  observeEvent(input$display_fusion, {
    controlsAll$display_fusion <- input$display_fusion
    controlsAll$refresh_table <- gv$FUSION_FLAG
    set_all_dfs() 
  })
  observeEvent(input$display_segmentation, {
    controlsAll$display_segmentation <- input$display_segmentation
    controlsAll$refresh_table <- gv$SEGMENTATION_FLAG
    set_all_dfs() 
  })
  observeEvent(input$display_sample, {
    cnb_echo(gv,output,flog.info,"input: observeEvent(input$display_sample): ")
    controlsAll$display_sample <- input$display_sample
    selected = (controlsAll$reference_samples!=controlsAll$display_sample)
    updateCheckboxInput(session, inputId = "ui_reference_samples", value = selected)
    controlsAll$refresh_table <-  gv$WG_FLAG + gv$TARGETED_FLAG + gv$OFF_TARGET_FLAG + gv$BAF_FLAG +gv$FUSION_FLAG + gv$SEGMENTATION_FLAG + gv$ERROR_BARS_FLAG
    if(init_table$first_time_flag==TRUE) {  
      init_table$first_time_flag = FALSE
    }
    set_all_dfs() 
  })

  # -------------------------------------------------------------------
  observeEvent(input$plot1_click, {
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$plot1_click): = ",input$plot1_click))
    click <- input$plot1_click
    displayOptions2$ymin <- displayOptions1$ymin
    displayOptions2$ymax <- controlsAll$cn_max_val
    x <- click$x
    idx <- which.min((x-gv$cyto.anno$start_p) +  1e12 * as.integer((x-gv$cyto.anno$start_p)<0))
    chr <- gv$cyto.anno[idx,] 
    displayOptions2$chr_label <- as.character(chr$chr_label)
    displayOptions2$chr <- as.character(chr$chr)
    if(x > chr$end_p) {
      displayOptions2$arm <- "q"
      displayOptions2$xmin <- chr$start_q
      displayOptions2$xmax <- chr$end_q
    } else {
      displayOptions2$arm <- "p"
      displayOptions2$xmin <- chr$start_p
      displayOptions2$xmax <- chr$end_p
    }
    displayOptions2$df <- selectPointsInGenomicRange(df$df2,displayOptions2,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
  })
  # -------------------------------------------------------------------
  observeEvent(input$plot2_brush,{
    displayOptions3$brush <- input$plot2_brush
    displayOptions3$brush$ymin <- displayOptions1$ymin # Set it to zero to save hassle with annotations etc 
    displayOptions3$ymin <- displayOptions3$brush$ymin
    displayOptions3$ymax <- displayOptions3$brush$ymax
    displayOptions3$xmin <- displayOptions3$brush$xmin 
    displayOptions3$xmax <- displayOptions3$brush$xmax
    displayOptions3$df <- selectPointsInBox(df$df3,displayOptions3,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
    displayOptions3$gene_names <- selectGenesInGenomicRange(gv$gene_names,displayOptions3,controlsAll$display_all_gene_names)
    displayOptions3$exons <- selectGenesInGenomicRange(gv$exons,displayOptions3,controlsAll$display_all_gene_names)
    displayOptions3$is_gene <- FALSE
    displayOptions3$chr <- displayOptions2$chr
    displayOptions3$chr_label <- displayOptions2$chr_label
    #     cnb_echo(gv,output,flog.info,paste0("displayOptions1$chr = ",displayOptions1$chr))
#     cnb_echo(gv,output,flog.info,paste0("displayOptions2$chr = ",displayOptions2$chr))
#     cnb_echo(gv,output,flog.info,paste0("displayOptions3$chr = ",displayOptions3$chr))
  })
  # -------------------------------------------------------------------
  observeEvent(input$plot3_brush,{
    brush <- input$plot3_brush
    df <- selectPointsInBox(df$df3,brush,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,NULL)
    # cnb_echo(gv,output,flog.info("(x1,x2,y1,y2)= (%f,%f,%f,%f)",brush$xmin,brush$xmax,brush$ymin,brush$ymax))
    # cnb_echo(gv,output,flog.info("Selected table for display: %d x %d",nrow(df),ncol(df)))
    if(nrow(df)>0) {
      idx  <- colSums(!is.na(df))!=0 # remove uninformative columns
      displayOptions4$df <- df[,idx]
      idx <- bitwAnd(df$type,gv$COUNTS_MASK) # replace type field with human-readable text
      idx <- ifelse(idx==gv$WG_FLAG,gv$FINE_IDX,idx) # TODO: Don't use WG_FLAG to identify count type
      displayOptions4$df$resolution <- ifelse(idx==0,"",gv$resolutions[idx])
      idx <- bitwAnd(df$type,gv$FILE_TYPE_MASK)
      displayOptions4$df$type <-  ifelse(idx==0,"counts",gv$file_type[idx])
      fapply_nearest_gene <- function(x,y) { which.min((x-y)^2)  }
      idx <- sapply(displayOptions4$df$x,fapply_nearest_gene,gv$exons$x)
      displayOptions4$df$N <- round(displayOptions4$df$N,2)
      displayOptions4$df$raw <- round(displayOptions4$df$raw,2)
      displayOptions4$df$gene <-  gv$exons[idx,"name"]
      displayOptions4$df$exon <-  gv$exons[idx,"exon_number"]
      displayOptions4$df$strand <-  gv$exons[idx,"strand"]
      displayOptions4$df$x <- NULL
      # Order for ease-of-use
      colorder <- 1:ncol(displayOptions4$df)
      names(colorder) <- names(displayOptions4$df)
      end_pos <- colorder["end"] # offset columns from end of genomic range col
      colorder["gene"] <- end_pos + 0.1
      colorder["exon"] <- end_pos + 0.2
      colorder["strand"] <- end_pos + 0.3
      colorder["type"] <- end_pos + 0.4
      if(sum(names(colorder)=="resolution")==1) { colorder["resolution"] <- end_pos + 0.5 }
      colorder["N"] <- end_pos + 0.6
      if(sum(names(colorder)=="Nsd")==1) { 
        displayOptions4$df$Nsd <- round(displayOptions4$df$Nsd,2)
        colorder["Nsd"] <- end_pos + 0.7 
      }
      displayOptions4$df <- displayOptions4$df[order(displayOptions4$df$start),order(colorder)]
    } else {
      displayOptions4$df <- NULL
    }
  })
  # -------------------------------------------------------------------
  observeEvent(input$plot3_dblclick, {
    cnb_echo(gv,output,flog.info,paste0("observeEvent(input$plot3_dblclick): = ",input$plot3_dblclick))
    dblclick <- input$plot3_dblclick
    if(displayOptions3$is_gene == TRUE) {
      displayOptions3$is_gene <- FALSE
      displayOptions3$ymin <- displayOptions3$brush$ymin
      displayOptions3$ymax <- displayOptions3$brush$ymax
      displayOptions3$xmin <- displayOptions3$brush$xmin
      displayOptions3$xmax <- displayOptions3$brush$xmax
      displayOptions3$df <- selectPointsInBox(df$df3,displayOptions3,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
      displayOptions3$gene_names <- selectGenesInGenomicRange(gv$gene_names,displayOptions3,controlsAll$display_all_gene_names)
      displayOptions3$exons <- selectGenesInGenomicRange(gv$exons,displayOptions3,controlsAll$display_all_gene_names)
    } else {
      displayOptions3$is_gene <- TRUE
      if(controlsAll$display_all_gene_names==TRUE) {
        displayOptions3$gene_names <- gv$gene_names[which.min((gv$gene_names$x-dblclick$x)^2),]
      } else {
        displayOptions3$gene_names <- gv$targeted_gene_names[which.min((gv$targeted_gene_names$x-dblclick$x)^2),]
      }
      cnb_echo(gv,output,flog.info,paste0("Looking for ",displayOptions3$gene_names$name," on chr ",displayOptions3$gene_names$chr))
      displayOptions3$exons <- gv$exons[(gv$exons$name==displayOptions3$gene_names$name),]
      displayOptions3$xmin <- min(displayOptions3$exons$x - displayOptions3$exons$len)
      displayOptions3$xmax <- max(displayOptions3$exons$x + displayOptions3$exons$len)
      displayOptions3$df <- selectPointsInBox(df$df3,displayOptions3,controlsAll$display_data_bitflag,controlsAll$min_probe_coverage,gv$hidden_points)
    }
  })
  # -------------------------------------------------------------------

  cnb_echo(gv,output,flog.info,"Render Plots")

  # -------------------------------------------------------------------
  plot_segmentation <- function(ret,displayOptions){
    if(controlsAll$display_segmentation) {
      idx <- (displayOptions$df$N<displayOptions$ymax) & (bitwAnd(displayOptions$df$type,gv$SEGMENTATION_FLAG)>0)
      if(sum(idx)>0) {
        text_offset <- gv$text_offset_fraction * (displayOptions$ymax - displayOptions$ymin)
        df <- displayOptions$df[idx,c("x","x1","x2","N","Ns")]   
        df$x1 <- ifelse(df$x1 < displayOptions$xmin,displayOptions$xmin,df$x1)
        df$x2 <- ifelse(df$x2 > displayOptions$xmax,displayOptions$xmax,df$x2)
        if(nrow(df)<=24) {
          ret <- ret + annotate("text",x = df$x1,y = df$N + text_offset,label = as.character(round(df$N,2)),size = 4,colour = "black")
        }
        return(ret + annotate("segment",x = df$x1,xend = df$x2,y = df$N,yend = df$N,size = 0.6,linetype = 2, alpha = 1,colour = gv$segmentation_colour))
      }
    }
    return(ret)
  }
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  plot_fusion <- function(ret,displayOptions){
    if(controlsAll$display_fusion) {  
      idx <- (bitwAnd(displayOptions$df$type,gv$FUSION_FLAG)>0)
      if(sum(idx)>0) {
        text_offset <- gv$text_offset_fraction * (displayOptions$ymax - displayOptions$ymin)
        df <- displayOptions$df[idx,c("x","x1","x2","N","nSupport","chr","chr2")]   
        df$N <- ifelse(df$N < displayOptions$ymax, df$N, displayOptions$ymax)
        idx1 <- (df$x1 >= displayOptions$xmin)
        idx2 <- (df$x2 <= displayOptions$xmax)
        if(sum(idx1)>0) { ret <- ret + annotate("point",x = df$x1[idx1],y = df$N[idx1],shape = 1,size = 4,colour = gv$fusion_colour ) }
        if(sum(idx2)>0) { ret <- ret + annotate("point",x = df$x2[idx2],y = df$N[idx2],shape = 3,size = 4,colour = gv$fusion_colour ) }
        df$x1 <- ifelse(df$x1 < displayOptions$xmin,displayOptions$xmin,df$x1)
        df$x2 <- ifelse(df$x2 > displayOptions$xmax,displayOptions$xmax,df$x2)
        # idx <- ((df$chr!=df$chr2) & (idx1 | idx2))
        idx <- (idx1 | idx2)
        if(sum(idx)>0) {
          txt <- gsub("chr","",paste0("t(",df$chr[idx],":",df$chr2[idx],") N = ",df$nSupport[idx]))
          ret <- ret + annotate("text",x = df$x1[idx],y = df$N[idx] + text_offset,label = txt,size = 4,colour = gv$fusion_colour )
          ret <- ret + annotate("segment",x = df$x1[idx],xend = df$x2[idx],y = df$N[idx],yend = df$N[idx],linetype = 1,size = gv$fusion_line_size,alpha = 1,colour = gv$fusion_colour)
        }
      }
    }
    return(ret)
  }
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  plot_blank <- function(displayOptions){
    df <- data.frame(Ns = 0:(gv$cn_different_colours-1),N = 0,x = 0,y=0,type=0)
    df$N[1:4] <- c(displayOptions$ymin,displayOptions$ymax,0,1)
    df$x[1:4] <- c(displayOptions$xmin,displayOptions$xmax,displayOptions$xmin,displayOptions$xmax)
    df$type[1:4] <- c(1,1,gv$BAF_FLAG,gv$BAF_FLAG)
    df$CN <- as.factor(df$Ns)
    return(geom_point(data=df,mapping=aes(x = x, y = N, colour = CN),alpha = 0))
  }
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  # Render Plots  
  # -------------------------------------------------------------------
  output$plot1 <- renderPlot({
    cnb_echo(gv,output,flog.info,"Render Plot1")
    if(is.null(displayOptions1$df)){return(qplot(x,y, data = data.frame(x=double(0),y=double(0)), geom = "blank") + xlab("") + ylab(""))}
    mask <- bitwOr(gv$COUNTS_MASK,gv$BAF_FLAG)
    idx <- (displayOptions1$df$N < displayOptions1$ymax) & (bitwAnd(displayOptions1$df$type,mask)>0)
    df_counts <- displayOptions1$df[idx,c("x","N","Ns")]
    df_counts$CN <- as.factor(df_counts$Ns)
    ret <- ggplot(df_counts, aes(x = x, y = N,size=0, colour = CN )) +
      # facet_grid(type ~ .) + 
      scale_colour_manual(values = gv$cn_colour_scale, guide="legend") +
      ggtitle(controlsAll$display_sample) + 
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank()
      ) + 
      scale_size_identity(aes(size=0), guide="none") + 
      ylim(displayOptions1$ymin,displayOptions1$ymax) +
      xlab("") + 
      ylab("N") +
      annotate("text", x = gv$cyto.anno$x, y = rep(controlsAll$cn_max_val,nrow(gv$cyto.anno)), label = gv$cyto.anno$chr_label)
    if(nrow(df_counts)>0) { ret <- ret +  geom_point() }
    ret <- plot_segmentation(ret,displayOptions1)
    ret <- plot_fusion(ret,displayOptions1)
    # ret <- ret + plot_blank(displayOptions1)
    df_cyto <- gv$cyto.anno
    df_cyto$CN <- as.factor(df_cyto$Ns)
    ret  <- ret +
      geom_rect(data=df_cyto,mapping=aes(xmin = start_p, xmax = end_p, ymin = -Inf, ymax = Inf), alpha = 0.1,fill = "blue", show.legend=FALSE) +
      geom_rect(data=df_cyto,mapping=aes(xmin = start_q, xmax = end_q, ymin = -Inf, ymax = Inf), alpha = 0.1,fill = "red", show.legend=FALSE)
    idx <- (displayOptions1$df$N > displayOptions1$ymax) & (bitwAnd(displayOptions1$df$type,gv$COUNTS_MASK)>0)
    if(sum(as.integer(idx))>0) {
      df_overflow <- displayOptions1$df[idx,c("x","Ns","N")]
      df_overflow$Nlim <- displayOptions1$ymax
      df_overflow$CN <- as.factor(df_overflow$Ns)
      ret <- ret + annotate("point",x = df_overflow$x,y = df_overflow$Nlim,size = gv$init_tab$plot1_point_size,colour = "red",shape = 2)
    }
    ret  
  })
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  output$plot2 <- renderPlot({
    cnb_echo(gv,output,flog.info,"Render Plot2")
    if(is.null(displayOptions2$df)) { return(qplot(x,y, data = data.frame(x=double(0),y=double(0)), geom = "blank") + xlab("") + ylab(""))}
    mask <- bitwOr(gv$COUNTS_MASK,gv$BAF_FLAG)
    idx <- (displayOptions2$df$N<displayOptions2$ymax) & (bitwAnd(displayOptions2$df$type,mask)>0)

    df_counts <- displayOptions2$df[idx,c("x","N","Ns")]
    df_counts$CN <- as.factor(df_counts$Ns)
    ret <- ggplot(df_counts, aes(x = x, y = N,size=gv$init_tab$plot2_point_size, colour = CN, xlab="", ylab="N")) +
      ggtitle(paste0(controlsAll$display_sample,": ",displayOptions2$chr_label, displayOptions2$arm)) + 
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        panel.background = element_rect(fill = gv$plot2_colour)
      ) + 
      scale_colour_manual(values = gv$cn_colour_scale, guide="legend") +
      scale_size_identity(aes(size=0), guide="none") + 
      xlim(displayOptions2$xmin,displayOptions2$xmax) +
      ylim(displayOptions2$ymin,displayOptions2$ymax) +
      xlab("") + 
      ylab("N")
    
    if(sum(as.integer(idx))>0) { ret <- ret +  geom_point() }
    idx <- (displayOptions2$df$N > displayOptions2$ymax) & (bitwAnd(displayOptions2$df$type,gv$COUNTS_MASK)>0)
    if(sum(as.integer(idx))>0) {
      overflow <- displayOptions2$df[idx,c("x","Ns")]
      overflow$Nlim <- displayOptions2$ymax
      overflow$text <- as.character(overflow$Ns)
      ret <- ret + annotate("point", x = overflow$x, y = overflow$Nlim, size = 4, colour = "red", shape = 2) 
      ret <- ret + annotate("text", x = overflow$x, y = overflow$Nlim + runif(nrow(overflow),min = -1,max = 0),label = overflow$text)
    }
    ret <- plot_segmentation(ret,displayOptions2)
    ret <- plot_fusion(ret,displayOptions2)
    # ret <- ret + plot_blank(displayOptions2)
    if(controlsAll$display_cytobands) {
      idx <- (gv$cyto.text$x1 >= displayOptions2$xmin) & (gv$cyto.text$x2 < displayOptions2$xmax)
      df <- gv$cyto.text[idx,]
      if(!is.null(df)) {
        df$N <- displayOptions2$ymax
        ret <- ret + 
          annotate(
            "segment",x = df$x1,xend = df$x2,y = df$N * 0.9,yend = df$N * 0.9,size = 3,colour = df$colour) +
            annotate("text",x = df$x1 +  df$len /2,y = df$N, label = df$name,size = 4,colour = "black"
          ) 
      }
    }
    return(ret)
  })
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  output$plot3 <- renderPlot({
    cnb_echo(gv,output,flog.info,"Render Plot3")
    if(is.null(displayOptions3$df)) {return(qplot(x,y, data = data.frame(x=double(0),y=double(0)), geom = "blank") + xlab("") + ylab(""))}
    offset <- gv$chrom.info[displayOptions3$chr,"total"]
    xmin <- round(displayOptions3$xmin-offset)
    xmax <- round(displayOptions3$xmax-offset)
    mask <- bitwOr(gv$COUNTS_MASK,gv$BAF_FLAG)
    idx <-  (bitwAnd(displayOptions3$df$type,mask)>0)
    if(sum(idx)==0) {return(qplot(x,y, data = data.frame(x=double(0),y=double(0)), geom = "blank") + xlab("") + ylab(""))}
    df_counts <- displayOptions3$df[idx,c("x","N","Ns","Nmin","Nmax")]
    df_counts$CN <- as.factor(df_counts$Ns)
    ret <- ggplot(df_counts, aes(x = x, y = N, size=gv$init_tab$plot3_point_size, colour = CN, xlab="", ylab="N")) +
      ggtitle(paste0(controlsAll$display_sample,"::",displayOptions3$chr_label,":",xmin,"-",xmax)) + 
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        # panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = gv$plot3_colour)
      ) + 
      scale_colour_manual(values = gv$cn_colour_scale, guide="legend") +
      scale_size_identity(aes(size=0), guide="none") + 
      xlim(displayOptions3$xmin,displayOptions3$xmax) +
      # ylim(as.character(round(displayOptions3$ymin):round(displayOptions3$ymax))) +
      ylim(displayOptions3$ymin,displayOptions3$ymax) +
      xlab("") + 
      ylab("N")
    if(nrow(df_counts)>0) {
      if(controlsAll$display_error_bars) { ret <- ret + geom_pointrange(aes(ymin = Nmin, ymax = Nmax)) } 
      else { ret <- ret + geom_point() }
    }    
    ret <- plot_segmentation(ret,displayOptions3)
    ret <- plot_fusion(ret,displayOptions3)
    # ret <- ret + plot_blank(displayOptions3)
    if(!is.null(displayOptions3$gene_names) &&  (nrow(displayOptions3$gene_names)>0)) {
      ret <- ret + annotate(
        "text", 
        x = displayOptions3$gene_names$x,
        y = min(displayOptions3$gene_names$N + runif(nrow(displayOptions3$gene_names),min = 0.4,max = 1),displayOptions3$ymax), 
        label = displayOptions3$gene_names$name
        )         
    }

    if(!is.null(displayOptions3$exons) &&  (nrow(displayOptions3$exons)>0) && !is.null(displayOptions3$exons$x)) {
      ret <- ret + 
      annotate("segment",x = displayOptions3$exons$x1,xend = displayOptions3$exons$x2,y = displayOptions3$exons$N,yend = displayOptions3$exons$N,size = 2,colour = "black") +
      annotate("point", x = displayOptions3$exons$x, y = displayOptions3$exons$N, size = 1, colour = "black") +
      annotate("text",x = displayOptions3$exons$x, y =  displayOptions3$exons$N + 0.2,label = displayOptions3$exons$exon_number)
    }
    if(gv$init_tab$debug || gv$init_tab$verbose) {  output$error_message_text <- renderText({ gv$log }) }
    return(ret)
  })
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  # Render text
  # -------------------------------------------------------------------
  # Generate a summary of the dataset
  progress$inc(amount = 0.1, message = "Rendering Text", detail = NULL)
  output$hilighted_data_summary <- renderPrint(summary(displayOptions4$df$N,digits=3))
  output$hilighted_data <- renderDataTable(displayOptions4$df,options = list(pageLength = 100))
  progress$close()
})

