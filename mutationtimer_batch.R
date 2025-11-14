## ========= 用户需改的 3 处路径 =========
DIR_SNV_VCF   <- "~/Desktop/RA/PDAC/Result/somatic_vcf"            # SNV/MNV/ID VCF.gz
DIR_CNV_TSV   <- "~/Desktop/RA/PDAC/Code/WGS/evolution/PhylogicNDT/phylogic_inputs/cnv_with_purity"                 # Each sample's CNV+purity TSV directory
DIR_SV_TABLE  <- "/Users/jiajunzhan/Desktop/RA/PDAC/Result/manta_sv"                        
OUT_DIR       <- "~/Desktop/RA/PDAC/Result/mutationtimer"                                    
META_TSV <- "~/Desktop/RA/PDAC/Result/summary_with_metadata.tsv"
# Sample name determination logic: Extracted from the VCF file name, for example, L2500897.pass.vcf.gz -> L2500897
SAMPLE_FROM_VCF <- function(f) sub("-somatic-PASS\\.vcf\\.gz$", "", basename(f))
## ========= Dependency package =========
suppressPackageStartupMessages({
  library(MutationTimeR)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(RColorBrewer)
  library(dplyr)
})


meta <- read.delim(META_TSV, stringsAsFactors = FALSE)

## ========= Extract the function: gender + WGD from master =========
get_gender_wgd <- function(sm, meta, bb) {
  # Match by sample_id; if your master is not this column name, change it here
  row <- meta[meta$Sample == sm, , drop = FALSE]
  if (nrow(row) > 0 && "Gender" %in% colnames(row)) {
    sx <- as.character(row$Gender[1])
    if (sx == "MALE")        Gender <- "male"
    if (sx == "FEMALE")    Gender <- "female"
  }
  if (nrow(row) > 0 && "WGD" %in% colnames(row)) {
    ws <- as.character(row$WGD[1])
    if (ws == "TRUE")  WGD <- TRUE
    if (ws == "FALSE") WGD <- FALSE
  }
  
  list(Gender = Gender, is_wgd = WGD)
}

# Analyze the breakend ALT direction, broadly categorized as h2h / t2t.
.classify_bnd_dir <- function(alt_one){
  if (is.na(alt_one)) return(NA_character_)
  a <- as.character(alt_one)
  if (grepl("^\\].*\\]$", a) || grepl("^\\[.*\\[$", a)) return("h2hINV")
  if (grepl("^\\].*\\[$", a) || grepl("^\\[.*\\]$", a)) return("t2tINV")
  NA_character_
}

# Convert the "long table TSV" into a .plotSv usable VCF object (only needs chrom/start/end/svtype/ID/MATEID/ALT).
table_to_sv_for_plot <- function(df, seqinfo_ref=NULL, style=c("NCBI","UCSC")){
  style <- match.arg(style)
  need_cols <- c("chrom","start","end","svtype","ID","MATEID","ALT")
  stopifnot(all(need_cols %in% names(df)))
  
  rr <- GRanges(df$chrom, IRanges(as.integer(df$start), width=1))
  suppressWarnings(seqlevelsStyle(rr) <- style)
  
  SVTYPE    <- as.character(df$svtype)
  MATECHROM <- rep(NA_character_, nrow(df))
  MATEPOS   <- rep(NA_integer_,  nrow(df))
  SVCLASS   <- rep(NA_character_, nrow(df))
  
  # DEL / DUP: Mate in the same chorm, position = end
  w_end <- which(SVTYPE %in% c("DEL","DUP") & !is.na(df$end))
  if (length(w_end)) {
    MATECHROM[w_end] <- as.character(df$chrom[w_end])
    MATEPOS[w_end]   <- as.integer(df$end[w_end])
    SVCLASS[w_end]   <- SVTYPE[w_end]
  }
  # INV (single record with END)
  w_inv <- which(SVTYPE == "INV" & !is.na(df$end))
  if (length(w_inv)) {
    MATECHROM[w_inv] <- as.character(df$chrom[w_inv])
    MATEPOS[w_inv]   <- as.integer(df$end[w_inv])
    alt_txt <- as.character(df$ALT[w_inv])
    sub_cls <- vapply(alt_txt, .classify_bnd_dir, character(1))
    sub_cls[is.na(sub_cls)] <- "h2hINV"
    SVCLASS[w_inv] <- sub_cls
  }
  # BND / SGL / TRX / INTRX：Use MATEID to find a partner
  id2idx  <- seq_len(nrow(df)); names(id2idx) <- as.character(df$ID)
  mate_id <- as.character(df$MATEID)
  w_bnd <- which(SVTYPE %in% c("BND","SGL","TRX","INTRX") & !is.na(mate_id) & mate_id %in% names(id2idx))
  if (length(w_bnd)) {
    mate_idx <- id2idx[mate_id[w_bnd]]
    MATECHROM[w_bnd] <- as.character(df$chrom[mate_idx])
    MATEPOS[w_bnd]   <- as.integer(df$start[mate_idx])
    SVCLASS[w_bnd]   <- "TRA"  # Default TRA
    
    same_chr <- w_bnd[ df$chrom[w_bnd] == df$chrom[mate_idx] ]
    if (length(same_chr)) {
      alt_str <- as.character(df$ALT[same_chr])
      inv_sub <- vapply(alt_str, .classify_bnd_dir, character(1))
      idx_h2h <- same_chr[which(inv_sub == "h2hINV")]
      idx_t2t <- same_chr[which(inv_sub == "t2tINV")]
      if (length(idx_h2h)) SVCLASS[idx_h2h] <- "h2hINV"
      if (length(idx_t2t)) SVCLASS[idx_t2t] <- "t2tINV"
    }
  }
  
  infoDF <- S4Vectors::DataFrame(
    SVTYPE    = SVTYPE,
    MATECHROM = MATECHROM,
    MATEPOS   = MATEPOS,
    SVCLASS   = SVCLASS
  )
  ok <- !is.na(infoDF$MATECHROM) & !is.na(infoDF$MATEPOS) & !is.na(infoDF$SVCLASS)
  if (!any(ok)) {
    message("  [SV] No bow-shaped figure to draw (mate coordinates or type incomplete)")
    v <- VCF(rowRanges=rr[FALSE], colData=S4Vectors::DataFrame(row.names="SAMPLE"))
    info(v) <- infoDF[FALSE,]
    info(v)$SVCLASS <- factor(character(0), levels=c("DEL","DUP","h2hINV","t2tINV","TRA"))
    return(v)
  }
  rr_ok   <- rr[ok]
  info_ok <- infoDF[ok, , drop=FALSE]
  
  v <- VCF(rowRanges=rr_ok, colData=S4Vectors::DataFrame(row.names="SAMPLE"))
  info(v) <- info_ok
  info(v)$SVCLASS <- factor(info(v)$SVCLASS, levels=c("DEL","DUP","h2hINV","t2tINV","TRA"))
  
  if (!is.null(seqinfo_ref)) {
    suppressWarnings(seqlevelsStyle(seqinfo_ref) <- style)
    keep <- intersect(seqlevels(v), seqlevels(seqinfo_ref))
    v <- keepSeqlevels(v, keep, pruning.mode="coarse")
    genome(v) <- genome(seqinfo_ref)
  }
  v
}

# Construct my_regions, refLengths, chrOffset (for use in plotSample)
make_ref_and_regions <- function(si_ref, style=c("NCBI","UCSC")){
  style <- match.arg(style)
  suppressWarnings(seqlevelsStyle(si_ref) <- style)
  # Only use 1:22, X, Y (if available)
  wanted <- intersect(seqlevels(si_ref), c(as.character(1:22),"X"))
  si_ref <- si_ref[wanted]
  refLengths <- GRanges(seqnames=names(seqlengths(si_ref)),
                        ranges=IRanges(1, seqlengths(si_ref)),
                        seqinfo=si_ref)
  # chrOffset: In the order of NCBI style (1..22, X, Y, MT)
  ord <- intersect(c(as.character(1:22),"X","Y"), names(seqlengths(si_ref)))
  lvec <- as.numeric(seqlengths(si_ref)[ord])
  chrOffset <- cumsum(c(0, lvec))
  names(chrOffset) <- c(ord, "MT")
  list(refLengths=refLengths, chrOffset=chrOffset,
       my_regions=as(si_ref[ord], "GRanges"))
}

# Save a sample drawing (PDF or PNG)
save_plot_sample <- function(file, expr, width=12, height=16, type=c("pdf","png")){
  type <- match.arg(type)
  if (type=="pdf") {
    pdf(file, width=width, height=height)
    on.exit(dev.off(), add=TRUE)
    eval.parent(substitute(expr))
  } else {
    png(file, width=width, height=height, units="in", res=180)
    on.exit(dev.off(), add=TRUE)
    eval.parent(substitute(expr))
  }
}

## ========= Create output directory =========
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "rdata"), showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(OUT_DIR, "batch.log")

## ========= List all VCF (as a sample list) =========
vcf_files <- Sys.glob(file.path(DIR_SNV_VCF, "*-somatic-PASS.vcf.gz"))
stopifnot(length(vcf_files) > 0)

## ========= Batch main loop =========
sink(log_file, split=TRUE)
cat(sprintf("[%s] Batch start: found %d VCF files\n", Sys.time(), length(vcf_files)))

for (f in vcf_files) {
  cat("\n", strrep("=", 80), "\n", sep="")
  sm <- SAMPLE_FROM_VCF(f)
  cat(sprintf("[%s] Sample %s\n", Sys.time(), sm))
  out_pdf <- file.path(OUT_DIR, "plots", paste0(sm, ".timing.pdf"))
  out_rds <- file.path(OUT_DIR, "rdata", paste0(sm, ".objects.rds"))
  # If both PDF and RDS exist, skip.
  if (file.exists(out_pdf) && file.exists(out_rds)) {
    cat(sprintf("[%s] Sample %s is done. Skip\n", Sys.time(), sm))
    next
  }
  
  try({
    # ---- 1) Read VCF, take tumor sample column ----
    vcf <- readVcf(f)
    vcf <- expand(vcf)
    # The default takes the second column as the tumor; if you need to specify manually, change it here:
    tumor_sample <- colnames(vcf)[2]
    vcf <- vcf[, tumor_sample, drop=FALSE]
    
    # Unified to NCBI style, and only retain the main chormosomes
    suppressWarnings(seqlevelsStyle(vcf) <- "NCBI")
    keep_v <- intersect(seqlevels(vcf), as.character(c(1:22,"X")))
    vcf <- keepSeqlevels(vcf, keep_v, pruning.mode="coarse")
    
    # ---- 2) Read CNV+purity table, build bb ----
    cnf <- file.path(DIR_CNV_TSV, paste0(sm, "_with_purity.tsv"))
    if (!file.exists(cnf)) stop("CNV+purity File missing: ", cnf)
    seg <- read.delim(cnf, stringsAsFactors = FALSE)
    seg$chromosome <- sub("^chr", "", seg$chromosome)
    bb <- GRanges(
      seqnames = seg$chromosome,
      ranges   = IRanges(seg$start, seg$end),
      major_cn = seg$majorAlleleCopyNumber,
      minor_cn = seg$minorAlleleCopyNumber,
      clonal_frequency = seg$purity
    )
    maj <- as.numeric(mcols(bb)$major_cn)
    min <- as.numeric(mcols(bb)$minor_cn)
    maj[maj < 0] <- 0; min[min < 0] <- 0
    maj2 <- maj
    min2 <- min
    mcols(bb)$major_cn <- maj2
    mcols(bb)$minor_cn <- min2
    mcols(bb)$major_cn <- as.integer(round(mcols(bb)$major_cn))
    mcols(bb)$minor_cn <- as.integer(round(mcols(bb)$minor_cn))
    
    # Align seqinfo
    si_ref <- seqinfo(rowRanges(vcf))
    suppressWarnings(seqlevelsStyle(bb) <- "NCBI")
    keep_b <- intersect(seqlevels(bb), seqlevels(si_ref))
    bb  <- keepSeqlevels(bb, keep_b, pruning.mode="coarse")
    genome(bb) <- genome(si_ref)
    if (all(is.na(seqlengths(bb)))) seqlengths(bb) <- seqlengths(si_ref)
    
    # ---- 3) MutationTimeR ----
    gw <- get_gender_wgd(sm, meta, bb)
    gender_now <- gw$Gender
    is_wgd_now <- gw$is_wgd
    
    cat(sprintf("  [META] Gender=%s, is_wgd=%s\n",
                gender_now,
                ifelse(is_wgd_now, "TRUE", "FALSE")))
 #   clusters <- data.frame(
 #     cluster    = 1,
 #     n_ssms     = nrow(vcf),
 #     proportion = unique(seg$purity)   # 或者 max(cn$clonal_frequency, na.rm=TRUE)
#    )
    mt <- mutationTime(
      vcf,
      bb,
      n.boot = 50,
      isWgd  = is_wgd_now,
      gender = gender_now, clusters=clusters
    )
    
    vcf <- addMutTime(vcf, mt$V)
    mcols(bb) <- cbind(mcols(bb), mt$T)
    mcols(bb)$type <- as.character(mcols(bb)$type)
    
    # Unified category level (for plotSample)
    info(vcf)$CLS <- factor(as.character(info(vcf)$CLS),
                            levels=c("clonal [early]","clonal [late]","clonal [NA]","subclonal"))
    bb$type <- factor(as.character(mcols(bb)$type),
                      levels=c("Mono-allelic Gain","CN-LOH","Bi-allelic Gain (WGD)"))
    
    # ---- 4) Read SV table and convert it into drawable objects ----
    sv_tab <- Sys.glob(file.path(DIR_SV_TABLE, paste0("*", sm, "*.tsv")))
    if (length(sv_tab) > 0) {
      df <- read.delim(sv_tab[1], stringsAsFactors = FALSE)
      if ("chrom" %in% names(df)) df$chrom <- sub("^chr", "", df$chrom)
      df <- df %>%
        mutate(tier = as.integer(trimws(as.character(tier)))) %>%
        filter(tier %in% c(1, 2, 3))
      sv_for_plot <- table_to_sv_for_plot(df, seqinfo_ref = si_ref, style = "NCBI")
      cat(sprintf("  [SV] Total drawable: %d\n", nrow(sv_for_plot)))
    } else {
      cat("  [SV] No corresponding SV table found, skip SV drawing\n")
      sv_for_plot <- NULL
    }
    
    # ---- 5) my_regions / refLengths / chrOffset ----
    refobj <- make_ref_and_regions(seqinfo(rowRanges(vcf)), style="NCBI")
    refLengths <<- refobj$refLengths
    chrOffset  <<- refobj$chrOffset
    my_regions <- refobj$my_regions
    
    # ---- 6) Plot & save ----
    save_plot_sample(out_pdf, {
      if (is.null(sv_for_plot) || nrow(sv_for_plot)==0) {
        plotSample(vcf, bb, regions = my_regions)  # no SV
      } else {
        plotSample(vcf, bb, sv = sv_for_plot, regions = my_regions)
      }
    }, type="pdf",
    width = 7, height = 7)
    cat("  [PLOT] save：", out_pdf, "\n")
    
    # ---- 7) save obj----
    saveRDS(list(vcf=vcf, bb=bb,
                 sv=sv_for_plot,
                 mt_V=mt$V, mt_T=mt$T,
                 my_regions=my_regions),
            out_rds)
    cat("  [RDS]  save：", out_rds, "\n")
    
  }, silent = TRUE) -> res
  
  if (inherits(res, "try-error")) {
    cat(sprintf("  [ERROR] %s\n", as.character(res)))
  } else {
    cat(sprintf("[%s] Sample %s is done\n", Sys.time(), sm))
  }
}

cat(sprintf("\n[%s] All completed\n", Sys.time()))
sink()