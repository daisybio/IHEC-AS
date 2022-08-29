aggregate_hits <- function(merge_dt, name, hits, signal_gr, peak = FALSE){
  if (peak){
    # merge_dt[, paste(region_name, 'presence', sep = '_') := ID %in% from(hits)]
    tmp_dt <- data.table(ID=from(hits))[, .(peak_count=.N), by=ID]
    merge_dt[tmp_dt, on=.(ID), (paste(name, 'peak_count', sep = ';')):=peak_count]
    
    # region_gr is an already aggregated data table if peak is true and enhancers are aggregated
    # if (endsWith(region_name, 'enhancer')) {
    #   merge_dt[region_gr, on=.(ID), paste(region_name, 'percentage', sep = '_') := percentage]
    # } else {
    #   overlaps <- pintersect(region_gr[from(hits)], signal_gr[to(hits)])
    #   percentOverlap <- width(overlaps) / width(region_gr[from(hits)])
    #   merge_dt[data.table(ID = from(hits), percentage = percentOverlap)[, .(percentage = sum(percentage)), by = ID], on=.(ID), paste(region_name, 'percentage', sep = '_') := percentage]
    # }
    # merge_dt[get(paste(region_name, 'presence', sep = '_')) == FALSE, paste(region_name, 'percentage', sep = '_') := 0]
    # invisible(merge_dt)
  } else {
    tmp_dt <- data.table(ID = hits@from, score = signal_gr[hits@to]$score)[, sapply(aggregation_functions, function(agg_method){agg_method(score)}, simplify = FALSE), by = ID]
    merge_dt[tmp_dt, on=.(ID), (paste(name, names(aggregation_functions), sep = ';')):=mget(names(aggregation_functions))]
    # res_dt <- rbindlist(sapply(aggregation_functions, function(fun) tmp_dt[, .(value=fun(score)), by=ID], simplify = FALSE), idcol='aggregation', use.names = TRUE)
    # res_dt[, aggregation:=as.factor(aggregation)]
  }
  invisible(merge_dt)
}

aggregate_with_flank <- function(merge_dt, name, signal_gr, event_gr, cCRE_dt, cCREs, cCRE_gr, upstream_gr, downstream_gr, flank_size, peak = FALSE, make_cCREs = FALSE) {
  
  region_grs <- list(
    # promoter_gr, # promoter regions
    upstream_other_region = upstream_gr, # the whole region upstream
    upstream_tile = flank(event_gr, flank_size), # 150 bp upstream
    event_name = event_gr, # alternative region
    downstream_tile = flank(event_gr, flank_size, start = FALSE), # 150 bp downstream
    downstream_other_region = downstream_gr # the whole region downstream
  )
  
  hits_list <- sapply(region_grs, findOverlaps, subject = signal_gr, simplify = FALSE)
  # names(hits_list) <- c('upstream_other_region', 'upstream_tile', 'event_name', 'downstream_tile', 'downstream_other_region')
  
  # hits_list <- c(hits_list, cCRE_hits_list)
  
  # enhancer regions
  # first find enhancers that are at most 5kb from the alternative region
  cCRE_hits <- findOverlaps(event_gr, cCRE_gr, maxgap = 5000)
  # get the hits of the enhancers that were near alternative regions
  cCRE_signal_hits <- findOverlaps(cCRE_gr[to(cCRE_hits)], signal_gr)
  # get the id of the alternative region
  from_cCRE <- from(cCRE_hits)[cCRE_signal_hits@from]
  # get the id in the signal file
  to_cCRE <- to(cCRE_signal_hits)
  # event_cCRE_hits <- data.table(
  #   from= from_cCRE,
  #   to=to_cCRE,
  #   cCRE=factor(cCREs[to(cCRE_hits)[cCRE_signal_hits@from], accession])
  # )
  # cCRE_dt_list <- split(event_cCRE_hits, by = 'cCRE', keep_by=FALSE)
  # cCRE_hits_list <- sapply(cCRE_dt_list, function(dt) dt[, Hits(from, to, nLnode = max(from), nRnode = max(to))], simplify = FALSE)
  
  # cCRE_acc2idx <- cCREs[, seq(.N)]
  # names(cCRE_acc2idx) <- cCREs[, accession]
  #TODO: cCRE rethinking
  
  from_split <- split(from_cCRE, cCREs[to(cCRE_hits)[cCRE_signal_hits@from], cCRE_type])
  to_split <- split(to_cCRE, cCREs[to(cCRE_hits)[cCRE_signal_hits@from], cCRE_type])
  hits_enhancer <- mapply(function(from, to) {Hits(from, to, max(from), max(to))}, from_split, to_split, SIMPLIFY = FALSE)
  # if (length(cCRE_signal_hits) == 0 || length(cCRE_hits) == 0) hits_enhancer <- Hits()
  # else hits_enhancer <- Hits(from_cCRE, to_cCRE, max(from_cCRE), max(to_cCRE))
  hits_list <- c(hits_enhancer, hits_list)
  # names(hits_list)[1] <- 'enhancer'
  
  names(hits_list) <- paste(name, names(hits_list), sep = ';')
  
  # if (peak) {
  #   overlap_enhancers <- pintersect(enhancer_gr[cCRE_hits@to][cCRE_signal_hits@from], signal_gr[cCRE_signal_hits@to])
  #   overlap_dt <- data.table(ID = from_cCRE, 
  #                            overlap_width = width(overlap_enhancers), 
  #                            enhancer_width = width(enhancer_gr[cCRE_hits@to][cCRE_signal_hits@from]))[
  #                              , .(overlap_width = sum(overlap_width), enhancer_width = sum(enhancer_width)), by=ID][
  #                                , percentage := overlap_width/enhancer_width]
  #   region_grs <- c(list(overlap_dt), region_grs)
  # } else region_grs <- c(list(enhancer_gr[cCRE_hits@to]), region_grs)
  # aggregate cCREs separately
  if (make_cCREs)
  aggregate_hits(merge_dt = cCRE_dt, name = name, hits = findOverlaps(cCRE_gr, signal_gr), signal_gr = signal_gr, peak = FALSE)
  
  invisible(mapply(aggregate_hits, name=names(hits_list), hits = hits_list, #region_gr = region_grs, 
         MoreArgs = list(merge_dt=merge_dt, signal_gr = signal_gr, peak = peak)))
  # res_dt <- rbindlist(mapply(aggregate_hits, hits = hits_list, #region_gr = region_grs, 
  #                            MoreArgs = list(signal_gr = signal_gr, peak = peak), SIMPLIFY = FALSE), idcol = 'region', use.names = TRUE)
  # res_dt[, region:=as.factor(region)]
  # res_dt
  cCRE_dt
}


### Main aggregation function ----
aggregate_multiple_samples <- function(samples_to_consider, event_dt, event_gr, upstream_gr, downstream_gr, make_cCREs = FALSE){
  
  # annotation <- import('suppa_analysis/gencode.v29.annotation.gtf')
  # used_genes <- annotation[annotation$gene_id %in% event_dt[, gene_id] & annotation$type == 'gene', ]
  # event_dt[data.table(gene_id = used_genes$gene_id, gene_start = start(used_genes), gene_end = end(used_genes)), on = 'gene_id', c('gene_start', 'gene_end') := .(gene_start, gene_end)]
  # 
  # event_dt[strand == '-', distance_TSS:=gene_end - end(event_gr[ID])]
  # event_dt[strand != '-', distance_TSS:=start(event_gr[ID]) - gene_start]
  # 
  # event_dt[strand == '-', distance_TES:=start(event_gr[ID]) - gene_start]
  # event_dt[strand != '-', distance_TES:=gene_end - end(event_gr[ID])]
  
  # # prepare promoters
  # promoters <- promoters(annotation[annotation$gene_id %in% event_dt[, gene_id] & annotation$type == 'gene', ], upstream = 1000, downstream = 0)
  # event_dt[data.table(gene_id = promoters$gene_id, seqnames = as.vector(seqnames(promoters)), start_promoter = start(promoters), end_promoter = end(promoters)), on = 'gene_id', c('start_promoter', 'end_promoter') := .(start_promoter, end_promoter)]
  # promoter_gr <- event_dt[, GRanges(seqnames, ranges = IRanges(ifelse(is.na(start_promoter), 0, start_promoter), ifelse(is.na(end_promoter), 0, end_promoter)), strand)]
  # 
  # # compute distances to promoters
  # upstream <- start(event_gr) - end(promoter_gr) 
  # downstream <- start(promoter_gr) - end(event_gr) 
  # distance_TSS <- pmax(upstream, downstream)
  # 
  # # compute distance to TES
  # tes <- flank(annotation[annotation$gene_id %in% event_dt[, gene_id] & annotation$type == 'gene', ], 10, start = FALSE)
  # event_dt[data.table(gene_id = tes$gene_id, seqnames = as.vector(seqnames(tes)), start_tes = start(tes), end_tes = end(tes)), on = 'gene_id', c('start_tes', 'end_tes') := .(start_tes, end_tes)]
  # tes_gr <- event_dt[, GRanges(seqnames, ranges = IRanges(ifelse(is.na(start_tes), 0, start_tes), ifelse(is.na(end_tes), 0, end_tes)), strand)]
  # upstream <- start(event_gr) - end(tes_gr) 
  # downstream <- start(tes_gr) - end(event_gr) 
  # distance_TES <- pmax(upstream, downstream)
  
  # add annotated enhancers and promoters
  cCREs <- data.table::fread('data/GRCh38-cCREs.bed')
  names(cCREs) <- c('seqnames', 'start', 'end', 'some_id', 'accession', 'cCRE_type')
  cCRE_gr <- cCREs[, GRanges(seqnames = seqnames, IRanges(start = start, end = end))]
  # mcols(enhancers) <- NULL
  cCRE_hits <- findOverlaps(event_gr, cCRE_gr, maxgap = 5000, ignore.strand=TRUE)
  # enhancers <- enhancers[unique(to(cCRE_hits))]
  # distances_enhancers <- distanceToNearest(event_gr, enhancers)
  # promoters_annotated <- rtracklayer::import('data/hg38_fair+new_CAGE_peaks_phase1and2.bed')
  # preceding <- precede(event_gr, promoters_annotated)
  # names(preceding) <- event_dt[, ID]
  # promoters_annotated <- promoters_annotated[preceding[!is.na(preceding)]]
  # region to cover when importing bigwigs to cover all regions of interest
  regions_to_cover <- c(cCRE_gr[unique(to(cCRE_hits))],
                        #promoter_gr, 
                        upstream_gr, 
                        event_gr, 
                        downstream_gr, 
                        promoters(event_gr, upstream = flank_size, downstream = flank_size))
  
  # ihec <- samples_to_consider[1]
  # profvis::profvis(
  #
  ihec_list <- pbmcapply::pbmclapply(samples_to_consider, function(ihec) 
    {
      
      merge_dt <- event_dt[, .(ID)]
      cCRE_dt <- data.table(ID = seq.int(nrow(cCREs)))
      # first gather gene expression
      # separator_gene_ids <- '_and_'
      # signal_psi_dt <- event_dt[, .(ID, PSI = get(ihec))] #, gene_id)]
      # signal_psi_dt <- signal_psi_dt[, .(gene_id = unlist(tstrsplit(gene_id, separator_gene_ids, fixed=TRUE))), by=.(ID, PSI)]
      # signal_psi_dt[unique(gene_quants[EpiRR_no_version == ihec, .(gene_id, gene_tpm)]), on = .(gene_id), gene_expression := gene_tpm]
      # signal_psi_dt <- signal_psi_dt[, .(gene_id=paste(gene_id, collapse=separator_gene_ids), gene_expression=sum(gene_expression, na.rm = TRUE)), by=.(ID, PSI)]
      
      # signal_psi_dt[, gene_id := NULL]
      # setkey(signal_psi_dt, ID)
      
      # mark_names <- c()
      # signal_gr_list <- list()
      # peak_list <- c()
      
      # add wgbs
      wgbs_files <- list.files(wgbs_data_dir, ihec)
      if (length(wgbs_files) == 2) {
        wgbs <- Reduce(c, lapply(c('pos', 'neg'), function(strand) {
          wgbs_file <-
            list.files(
              wgbs_data_dir,
              pattern = paste0('.*', ihec, '.*.gembs_', strand, '.bw$'),
              full.names = TRUE
            )
          wgbs_strand <- rtracklayer::import(wgbs_file, which = regions_to_cover)
          strand(wgbs_strand) <- ifelse(strand == 'pos', '+', '-')
          return(wgbs_strand)
        }))
        # mark_names <- c(mark_names, 'wgbs')
        # signal_gr_list <- c(signal_gr_list, wgbs)
        # peak_list <- c(peak_list, FALSE)
        
        aggregate_with_flank(merge_dt, 'wgbs', wgbs, event_gr, cCRE_dt, cCREs, cCRE_gr, upstream_gr, downstream_gr, flank_size, make_cCREs = make_cCREs)
      } else if (length(wgbs_files) > 2)
        warning(paste(ihec, 'had > 2 files:', paste(wgbs_files, collapse = ', ')))
      else if (length(wgbs_files) < 0) {
        warning(paste(ihec, 'had < 2 files'))
      }
      # add histone marks
      histone_samples <- fread(file.path(data_dir, 'ihec_metadata.csv'))
      histone_samples[, epirr_id_wo_version := tstrsplit(epirr_id, '.', fixed = TRUE)[1]]
      for (this_uuid in histone_samples[epirr_id_wo_version == ihec, uuid]) {
        file_ext <- c(peak = '\\.pval0\\.01\\.500K\\.bfilt\\.narrowPeak\\.gz$', 
                      signal = '\\.fc\\.signal\\.bigwig$')
        for (file_type in names(file_ext)) {
          histone_file <-
            list.files(
              chip_data_dir,
              pattern = paste0('.*', ihec, '.*', this_uuid, file_ext[file_type]),
              full.names = TRUE
            )
          if (length(histone_file) > 1)
            warning(paste(ihec, 'had > 1 files for uuid', this_uuid, ':', paste(histone_file, collapse = ', ')))
          else if (length(histone_file) == 0) {
            warning(paste(ihec, 'had 0 files for uuid', this_uuid))
          } else {
            hPTM <- rtracklayer::import(histone_file, which = regions_to_cover)
            
            # mark_names <- c(mark_names, histone_samples[uuid == this_uuid, antibody])
            # signal_gr_list <- c(signal_gr_list, hPTM)
            # peak_list <- c(peak_list, file_type == 'peak')
            
            aggregate_with_flank(merge_dt, histone_samples[uuid == this_uuid, antibody], hPTM, event_gr, cCRE_dt, cCREs, cCRE_gr, upstream_gr, downstream_gr, flank_size, peak = file_type == 'peak',  make_cCREs = make_cCREs)
          }
        }
      }
      
      # browser()
      # mark_list <-
      #   mapply(
      #     aggregate_with_flank,
      #     signal_gr_list,
      #     peak = peak_list,
      #     MoreArgs = list(
      #       event_gr = event_gr,
      #       cCREs = cCREs,
      #       cCRE_gr = cCRE_gr,
      #       upstream_gr = upstream_gr,
      #       downstream_gr = downstream_gr,
      #       flank_size = flank_size
      #     ),
      #     SIMPLIFY = FALSE
      #   )
      # names(mark_list) <- mark_names
      # res_dt <- rbindlist(mark_list, idcol = 'mark', use.names = TRUE)NULL
      # res_dt[, mark:=as.factor(mark)]
      
      # return(res_dt)
      fwrite(merge_dt, file.path(psi_input_dir, 'sample_dts', paste0(ihec, '-merge_dt.csv.gz')))
      if (make_cCREs)
      fwrite(cCRE_dt[apply(cCRE_dt, 1, function(row) sum(is.na(row)) < length(cCRE_dt) - 1)], file.path(psi_input_dir, 'sample_dts', paste0(ihec, '-cCRE_dt.csv.gz')))
      
      # merge_dt
      # list(merge_dt = merge_dt, cCRE_dt = cCRE_dt)
      TRUE
    })
  # 
  
  names(ihec_list) <- samples_to_consider
  # res_dt <- rbindlist(ihec_list, idcol = 'ihec', use.names = TRUE)
  # res_dt[, ihec:= as.factor(ihec)]
  # return(res_dt)
  ihec_list
  
}