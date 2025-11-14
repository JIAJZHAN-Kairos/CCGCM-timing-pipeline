## PCAWG-11 Mutational Signatures timing —— CCGCM Version

## Input files
## Libraries
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(plyr)
  library(scales)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(nnls)
  library(reshape2)
})

## Inputs
args <- commandArgs(TRUE)
id           <- args[1]
obj_rds      <- args[2]  
sbs_path     <- args[3]  # COSMIC_v3.4_SBS_GRCh38.txt
dbs_path     <- args[4]  # COSMIC_v3.4_DBS_GRCh38.txt
out_prefix <- args[5]

#obj_rds <- "/Users/jiajunzhan/Desktop/RA/PDAC/Result/mutationtimer/rdata/ACTN01020020.objects.rds"
#id <- "ACTN01020020"
#sbs_path <- "/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/COSMIC_v3.4_SBS_GRCh38.txt"
#dbs_path  <- "~/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/COSMIC_v3.4_DBS_GRCh38.txt"
#id_path   <- "~/Desktop/RA/PDAC/Code/WGS/evolution/signature/PCAWG11-Timing_and_Signatures-master/COSMIC_v3.4_ID_GRCh37.txt"
#out_prefix <- "/Users/jiajunzhan/Desktop/RA/PDAC/Code/WGS/evolution/signature/output/ACTN01020020/"
#n_boot <- 200
if (!dir.exists(out_prefix)) dir.create(out_prefix, recursive = TRUE)
setwd(out_prefix)


message("ID = ", id)
obj <- readRDS(obj_rds)
vcf <- obj$vcf
mt  <- as.data.frame(obj$mt_V) # ( “clonal [early]”, “clonal [late]”, “clonal”, “subclonal”）


get_alt_width <- function(vcf) {
  a <- alt(vcf)
  # Normal Condition: DNAStringSetList
  if (inherits(a, "DNAStringSetList")) {
    return(S4Vectors::elementNROWS(a))
  }
  # Sometimes, readVcfAsVRanges may directly return DNAStringSet.
  if (inherits(a, "DNAStringSet")) {
    return(rep(1, length(a)))
  }
  # If the extreme case is DNAString (single ALT)
  if (inherits(a, "DNAString")) {
    return(rep(1, length(vcf)))
  }
  stop("Unexpected ALT type: ", class(a))
}

refW <- Biostrings::width(VariantAnnotation::ref(vcf))
altW <- get_alt_width(vcf)
idx_snv <- refW == 1 & altW == 1
idx_id  <- (refW != altW)
vcf_snv_mnv <- vcf[idx_snv]
vcf_indel  <- vcf[idx_id]

cls <- as.character(mt$CLS)

t_df <- data.frame(
  seqnames = as.character(GenomicRanges::seqnames(vcf)),
  position = GenomicRanges::start(vcf),
  mut_type = ifelse(idx_id, "indel", "SNV"),   # Note: DBS is also set to "SNV" (consistent with the original script).
  timing   = cls,
  stringsAsFactors = FALSE
)


## Functions

# Determines if individual SNVs are part of an MNV
isMNV <- function(vcf) {
  d <- diff(start(vcf)) == 1 # are they next to each other?
  w <- c(FALSE, d) | c(d, FALSE) # need to add comparisons at the end, set to FALSE
  return(w)
}

# Formatting of SNVs -> pyrimidine context
formatSNVs <- function(vcf){
  vcf=vcf
  k=3
  ranges = resize(vcf, k, fix = "center") # resize to get 1 bp either side 
  seqlevelsStyle(ranges) = "UCSC"
  ref=BSgenome.Hsapiens.UCSC.hg38
  context = getSeq(ref, ranges) # get 3 bp contexts
  ct = as.data.frame(context)
  vcf$context = ct$x
  
  rev.context = reverseComplement(DNAStringSet(ct$x)) # get reverse complement of context as well
  vcf$rev.context = strsplit(toString(rev.context), ", ")[[1]]
  
  # formatting of bp change and contexts
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "A") | (ref(vcf) == "G" & alt(vcf) == "T")] = "C>A"
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "G") | (ref(vcf) == "G" & alt(vcf) == "C")] = "C>G"
  vcf$change[(ref(vcf) == "C" & alt(vcf) == "T") | (ref(vcf) == "G" & alt(vcf) == "A")] = "C>T"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "A") | (ref(vcf) == "A" & alt(vcf) == "T")] = "T>A"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "C") | (ref(vcf) == "A" & alt(vcf) == "G")] = "T>C"
  vcf$change[(ref(vcf) == "T" & alt(vcf) == "G") | (ref(vcf) == "A" & alt(vcf) == "C")] = "T>G"
  
  pyrimidine = c("C", "T")
  
  vcf$base_before = ifelse(ref(vcf) %in% pyrimidine, substring(vcf$context, 1, 1), substring(vcf$rev.context, 1, 1))
  vcf$base_after = ifelse(ref(vcf) %in% pyrimidine, substring(vcf$context, 3), substring(vcf$rev.context, 3))
  
  vcf$tnc = paste0(vcf$base_before,"[",vcf$change,"]",vcf$base_after)
  return(vcf)
}

# Format MNVs so that they match the signatures annotation
formatMNVs <- function(vcf, t_df){

  u <- reduce(granges(vcf[which(isMNV(vcf)),]))
  u <- u[width(u)==2]
  
  if(length(u)==0) return(u)
  
  seqlevelsStyle(u) = "UCSC"
  
  # Check timing is the same
  t = t_df
  t$start = t$position
  t$end = t$position
  t = subset(t, mut_type != "SV" & mut_type != "indel")
  t = makeGRangesFromDataFrame(t, keep.extra.columns = TRUE)
  seqlevelsStyle(t) = "UCSC"
  
  # Merge VCF with timing annotation
  timed_mnvs = mergeByOverlaps(t, u)
  
  # Remove single variants
  timed_mnvs$lengths = rep(rle(start(timed_mnvs$u))$lengths, rle(start(timed_mnvs$u))$lengths)
  timed_mnvs = subset(timed_mnvs, lengths>1)
  
  if(nrow(timed_mnvs)==0) return(timed_mnvs$u)
  
  times = lapply(split(timed_mnvs$timing, sort(rep(1:(nrow(timed_mnvs)/2), 2))), unique)
  times = !lapply(times, length) > 1
  
  u = reduce(timed_mnvs$u)[times]
  
  if(length(u)==0) return(u)

  seqlevelsStyle(vcf) = "UCSC"
  r <- as.character(ref(vcf))[vcf %over% u] # get reference base
  h <- subjectHits(findOverlaps(granges(vcf),u))
  a <- as.character(alt(vcf))[vcf %over% u] # get alternate base
  rr <- sapply(split(r, h), paste, collapse="")
  aa <- sapply(split(a, h), paste, collapse="")
  
  u$ref_mnv = rr
  u$alt_mnv = aa
  
  convert_all = c("AA","AG","CA","GA","GG","GT")
  
  # change purine bases, this is so they match the bases given in DBS signatures file
  c = which(u$ref_mnv %in% convert_all)
  u$ref_mnv[c] = reverseComplement(DNAStringSet(u$ref_mnv[c])) # if there are no purine bases it just won't add anything
  u$alt_mnv[c] = reverseComplement(DNAStringSet(u$alt_mnv[c]))
  
  # Format
  u$DBS = paste0(u$ref_mnv, ">", u$alt_mnv)
  
  # There are specific cases that need to be converted
  u$DBS[u$DBS=="AT>GG"] = "AT>CC"
  u$DBS[u$DBS=="AT>TC"] = "AT>GA"
  u$DBS[u$DBS=="AT>TG"] = "AT>CA"
  u$DBS[u$DBS=="CG>AA"] = "CG>TT"
  u$DBS[u$DBS=="CG>AC"] = "CG>GT"
  u$DBS[u$DBS=="CG>GA"] = "CG>TC"
  u$DBS[u$DBS=="GC>CT"] = "GC>AG"
  u$DBS[u$DBS=="GC>TG"] = "GC>CA"
  u$DBS[u$DBS=="GC>TT"] = "GC>AA"
  u$DBS[u$DBS=="TA>AC"] = "TA>GT"
  u$DBS[u$DBS=="TA>AG"] = "TA>CT"
  u$DBS[u$DBS=="TA>CC"] = "TA>GG"
  
  return(u)
}

seqDiff <- function(x,y,diff){
  x1 = unlist(strsplit(x, ""))
  y1 = unlist(strsplit(y, "")) 
  if (diff == "deleted"){
    return(x1[!1:length(x1) %in% pmatch(y1, x1)])
  } 
  if (diff == "inserted"){
    return(y1[!1:length(y1) %in% pmatch(x1, y1)])
  }
}


# Get indels that correspond to characteristic events from indel signatures ID1, ID2, ID13 and ID8
formatIndels <- function(vcf){
  
  indel_vcf = vcf
  
  if(length(indel_vcf)<=1) return(indel_vcf)
  
  seqlevelsStyle(indel_vcf) = "UCSC"
  
  # Start by computing what exactly is changed in the indel
  # NB - Some are both insertions and deletions e.g. AAA -> TT
  
  # Are all the reference bases in the alternate allele, in the right order?
  indel_vcf$deleted = !mapply(grepl, ref(indel_vcf), alt(indel_vcf))
  indel_vcf$inserted = !mapply(grepl, alt(indel_vcf), ref(indel_vcf))
  indel_vcf$del_base = mapply(seqDiff, ref(indel_vcf), alt(indel_vcf), diff="deleted", SIMPLIFY = FALSE)
  indel_vcf$ins_base = mapply(seqDiff, ref(indel_vcf), alt(indel_vcf), diff="inserted", SIMPLIFY = FALSE)
  
  # We are looking for 4 different types of indels, corresponding to ID1, ID2, ID8 and ID13
  
  # 1. ID1: 1bp insertion of T at regions of 5+ T's, or A at region of 5+ A's
  ins_1bp_T = indel_vcf[indel_vcf$ins_base=="T" | indel_vcf$ins_base=="A"]
  
  if (length(ins_1bp_T)>0){
    
    # Get 11bp surrounding, rle of T's/ A's, are they longer or equal to 5? 
    k=11
    ins_1bp_T = resize(ins_1bp_T, k, fix = "center") # get 5 bp either side of reference
    ref=BSgenome.Hsapiens.UCSC.hg38
    context_ins_T = as.data.frame(getSeq(ref, ins_1bp_T))
    context_ins_T$x = substring(context_ins_T$x, 2) # take the first base off, so that 5bp sequences must include centre base or be on the right hand side
    
    ins_1bp_T_5hp = c()
    for (n in 1:length(ins_1bp_T)){
      r = rle(strsplit(context_ins_T$x[n], "")[[1]]) # go through and count how many consecutive bases
      if (any(r$lengths >= 5)){ # are there any stretches >= 5?
        b = r$values[which(r$lengths>=5)] # get the 5+ repeated base
        ins_b = ins_1bp_T$ins_base[n] # and the inserted base
        if (ins_b %in% b){ # is inserted bp one of the repeated ones? In this case, this is fine to be either
          ins_1bp_T_5hp = c(ins_1bp_T_5hp, n) # if they are (should be A and A or T and T, then add this number to vector)
        }
      } 
    }
    
    ID1 = ins_1bp_T[ins_1bp_T_5hp]
    
    if (length(ID1)>0){ 
      ID1$change = NULL
      k=1
      ID1 = resize(ID1, k, fix="center") # resize back to centre single position
      ID1$classification = "ID1" 
    } # otherwise it is just the empty granges to be added
    
  } else {
    
    ID1 = ins_1bp_T # this will be empty, can just add to the others at the end
    
  }
  
  # 2. ID2: 1bp deletion of T at regions of 5+ T's, or A at regions of 5+ A's
  del_1bp_T = indel_vcf[indel_vcf$del_base=="T" | indel_vcf$del_base=="A"]
  
  if (length(del_1bp_T)>0){
    
    # Get 11bp surrounding, rle of T's/ A's, are they longer or equal to 6? 
    k=12
    del_1bp_T = resize(del_1bp_T, k, fix = "center") 
    ref=BSgenome.Hsapiens.UCSC.hg38
    context_del_T = as.data.frame(getSeq(ref, del_1bp_T))
    context_del_T$x = substring(context_del_T$x, 2)
    
    del_1bp_T_6hp = c()
    for (n in 1:length(del_1bp_T)){
      r = rle(strsplit(context_del_T$x[n], "")[[1]]) # go through and count how many consecutive bases
      if (any(r$lengths >= 6)){ # are there any stretches >= 6?
        b = r$values[which(r$lengths>=6)] # get the 6+ repeated base
        if (length(b)>1) { b = b[2] }
        del_b = del_1bp_T$del_base[n] # and the deleted base
        if (del_b %in% b){ # are they the same?
          del_1bp_T_6hp = c(del_1bp_T_6hp, n) # if they are (should be A and A or T and T, then add this number to vector)
        }
      } 
    }
    
    ID2 = del_1bp_T[del_1bp_T_6hp]
    
    if (length(ID2)>0) {
      ID2$change = NULL
      k = 1
      ID2 = resize(ID2, k, fix="center")
      ID2$classification = "ID2"
    }
    
  } else {
    ID2 = del_1bp_T
  }
  
  
  # 3. ID13: Get 1bp deletions of T, homopolymer length 2 (again, can either be T or A)
  del_1bp_T = indel_vcf[indel_vcf$del_base=="T" | indel_vcf$del_base=="A"]
  if (length(del_1bp_T)>0) {
    
    # ID13 is 1bp deletion of A or T, at either AA or TT
    k = 1
    del_1bp_T_3tnc = resize(del_1bp_T, k, fix = "start")
    k = 7
    del_1bp_T_3tnc = resize(del_1bp_T_3tnc, k, fix="center")
    context_del_1bp_T_3tnc = as.data.frame(getSeq(ref, del_1bp_T_3tnc))
    
    # I assume that the deleted base is always on the right in the reference allele
    del_1bp_T_2hp = c()
    for (n in 1:length(del_1bp_T)){
      
      # check is that base repeated either at the 4, 5 position, or the 5, 6 position, and not at 3 or 7
      tbp = strsplit(substring(context_del_1bp_T_3tnc$x[n], 4, 5), "")[[1]]
      fbp = strsplit(substring(context_del_1bp_T_3tnc$x[n], 5, 6), "")[[1]]
      tbp3 = substring(context_del_1bp_T_3tnc$x[n], 3, 3)
      fbp7 = substring(context_del_1bp_T_3tnc$x[n], 7, 7)
      del_b = del_1bp_T_3tnc$del_base[n]
      
      if ((all(tbp==del_b) && tbp3!=del_b) | (all(fbp==del_b) && fbp7!=del_b) ) {
        del_1bp_T_2hp = c(del_1bp_T_2hp, n)
      }
    }
      
    ID13 = del_1bp_T[del_1bp_T_2hp]
      
    if (length(ID13)>0){
      ID13$classification = "ID13"
      k=1
      ID13 = resize(ID13, k, fix="center")
    }
    
    
  } else { ID13 = del_1bp_T } # i.e. it is empty
  
  
  # 4. 5 bp + deletions, where the deleted segment is not repeated, and there is at least 1bp of microhomology
  bp5r = indel_vcf[which(lapply(indel_vcf$del_base, length)>=5)]
  
  if (length(bp5r)>0){
    bp5r$id = seq(length(bp5r)) # assign IDs to variants in the old file so can map back after changing sequence ranges
    bp5 = bp5r
    
    # is the deleted sequence repeated?
    l = length(bp5$del_base) # how long is deleted segment?
    bp5$l = l
    bp5 = bp5 + l # Add the length of the deletion either side
    
    ref = BSgenome.Hsapiens.UCSC.hg38
    context_bp5 = as.data.frame(getSeq(ref, bp5)) # get this wide range
    context_bp5$x = substring(context_bp5$x, 2)
    
    ms = c() # count how many repeats of the deleted segment per surrounding region
    for (m in 1:length(bp5)){
      counts = length(gregexpr(bp5$del_base[m], context_bp5$x[m]))
    ms = c(ms, counts)
    }
    
    bp5$context = context_bp5$x
    bp5 = bp5[ms==1] # the deleted segment is not repeated in the surrounding region
    
    
    if (length(bp5)>0){ # if there are any regions that do have repeats
      
      # Now look for microhomology in the deletions which are not repeated
      # The first deleted base has to be the same as the next remaining base
      bp5$nextbase = substring(bp5$context, (2*bp5$l)+1, (2*bp5$l)+1)
      mh = bp5[bp5$del_base[1][[1]][1]==bp5$nextbase]
      
      ID8 = bp5r[bp5r$id %in% mh$id]
      
      if (length(ID8)>0){
        ID8$classification = "ID8"
        ID8$del = NULL
        ID8$l = NULL
        ID8$context = NULL
        ID8$nextbase = NULL
        ID8$id = NULL
      } 	
    } else { ID8 = bp5r }
  } else { ID8 = bp5r}
  
  total_IDs = c(ID1, ID2, ID13, ID8)
  return(total_IDs)
  
} 




# Get timed multinomials of SNVs, MNVs and indels, removing MNVs from SNV calls
getSNVs_MNVs_IDs <- function(vcf_snv_mnv, vcf_indel, t_df){
  
  vcf = as(vcf_snv_mnv, "VRanges")
  
  # Get MNVs from the SNV vcf
  mnv_vcf = vcf[which(isMNV(vcf))]
  
  # Remove MNVs from SNVs
  snv_vcf = subsetByOverlaps(vcf, mnv_vcf, invert=TRUE)
  
  # Now get indels
  indel_vcf = as(vcf_indel, "VRanges")
  # Format each type of event to agree with classification in signatures
  
  # SNVs
  snv_vcf = formatSNVs(vcf=snv_vcf)
  snv = GRanges(seqnames = seqnames(snv_vcf), ranges = ranges(snv_vcf), class=snv_vcf$tnc) # Simplify VCF to just position and mutation class
  seqlevelsStyle(snv) = "UCSC"
  snv$type = "SNV"
  
  # MNVs
  mnv_vcf = formatMNVs(vcf, t_df) # Give it the unfiltered SNV VCF, as MNVs are identified within function
  seqlevelsStyle(mnv_vcf) = "UCSC"
  if (length(mnv_vcf)>0){ # Otherwise just leave mnv_vcf as empty GRanges, will be added to total and contribute nothing
    mnv_vcf$ref_mnv = NULL
    mnv_vcf$alt_mnv = NULL
    mnv_vcf$class = mnv_vcf$DBS
    mnv_vcf$DBS = NULL
    mnv_vcf$type = "SNV"
  }  

  # Indels
  indel = formatIndels(vcf=indel_vcf)
  if(length(indel)>1){
    seqlevelsStyle(indel) = "UCSC"
    indel = GRanges(seqnames = seqnames(indel), ranges = ranges(indel), class=indel$classification)
    rm(indel_vcf)
    indel$type = "indel"
  } else {
    indel = GRanges()
  }

  # Combine
  total_vcf = c(snv, mnv_vcf, indel)
  rm(vcf)
  
  # Get timing annotation
  t = t_df
  t$start = t$position
  t$end = t$position
  t = subset(t, mut_type != "SV")
  t = makeGRangesFromDataFrame(t, keep.extra.columns = TRUE)
  seqlevelsStyle(t) = "UCSC"
  
  # Merge VCF with timing annotation
  timed = mergeByOverlaps(t, total_vcf, type="start")
  timed = subset(timed, mut_type==type)
  
  # Double "clonal" section and add
  if ("clonal" %in% timed$timing){
    clonal = subset(timed, timing!="subclonal")
    clonal$timing = "clonal"
    timed = rbind(timed, clonal)
    # Also add total
    all = subset(timed, timing!="clonal")
    all$timing = "total"
    timed = rbind(timed, all)
  } else{
    timed = timed
    all = timed
    all$timing = "total"
    timed = rbind(timed, all)
  }
  
  # Split by timing
  timed = split(timed, timed$timing)
  
  # Make multinomial comprising all classes, per time frame
  timed_multi = lapply(timed, function(x){
    table(x$class)
  
  })
  
}

# Applies NNLS
nnls_sol <- function(m, s_active, no_active) {
  
  a = matrix(0, nrow = 1, ncol = no_active)
  for (i in 1:ncol(m))
    a[i,] = coef(nnls(s_active, m[,i]))
  a  
}

# Get signature weights from mutation matrix
getSignatureWeights <- function(mutationMatrix, active_comp, sig_comp, sigs){
  
  m = as.matrix(mutationMatrix)
  
  if (is.matrix(m) & is.numeric(m[1]) & is.matrix(active_comp) & is.numeric(active_comp[1]) & nrow(active_comp) == nrow(m)){
    
    no_active = ncol(active_comp)
    a = nnls_sol(m,active_comp,no_active)
    
    colnames(a) = sigs
    
    # add missing signatures set to 0 
    a_df = as.data.frame(a)
    add = subset(colnames(sig_comp), !colnames(sig_comp) %in% colnames(a_df))
    vals = rep(0,length(add))
    names(vals) = add
    vals = t(as.data.frame(vals))
    rownames(vals) = NULL
    a_final = cbind(a, vals)
    a_final = as.data.frame(a_final)
    a_final = a_final[,colnames(sig_comp)]
    a_final$sample = id
    
    return(a_final)
    
  } 
}


extractSigs <- function(events){

  # start with SNV signatures
  active_sigs = colnames(snv_comp)
  active_comp = as.matrix(snv_comp[,active_sigs])
  events_snv = lapply(events, function(x){
    x = x[names(x) %in% rownames(snv_comp)]
    if (length(names(x))!=96){
      missing = rownames(snv_comp)[!rownames(snv_comp) %in% names(x)]
      add = rep(0, length(missing))
      names(add) = missing
      x = c(x, add)
      x = x[rownames(snv_comp)]
    }
    x = x[rownames(snv_comp)]
    return(x)
  })
  weightsList = lapply(events_snv, getSignatureWeights, active_comp=active_comp, sig_comp=snv_comp, sigs=active_sigs)
  weights_snv = do.call("rbind", weightsList)
  
  # MNV signatures
    active_mnv = colnames(mnv_comp)
    active_mnv_comp =  as.matrix(mnv_comp[,active_mnv])
    events_mnv = lapply(events, function(x){
      x = x[names(x) %in% rownames(mnv_comp)]
      if (length(names(x))!=78){
        missing = rownames(mnv_comp)[!rownames(mnv_comp) %in% names(x)]
        add = rep(0, length(missing))
        names(add) = missing
        x = c(x, add)
        x = x[rownames(mnv_comp)]
      } 
      x = x[rownames(mnv_comp)]
      return(x)
    })
    weightsList_MNV = lapply(events_mnv, getSignatureWeights,active=active_mnv_comp,sig_comp=mnv_comp,sig=active_mnv)
    weights_mnv = do.call("rbind", weightsList_MNV)
  
  # ID signatures
  events_id = lapply(events, function(x){
    x = x[grep("ID", names(x))]
    indel_sigs = c("ID1", "ID2", "ID13", "ID8")
    missing = indel_sigs[which(!indel_sigs %in% names(x))]
    add = rep(0, length(missing))
    names(add) = missing
    x = c(x, add)
    x = x[indel_sigs]
  })
  events_id = do.call("rbind", events_id)
  if(ncol(events_id)!=4){
    add = data.frame(matrix(0, nrow=nrow(events_id), ncol=length(missing)))
    colnames(add) = missing
    events_id = cbind(events_id, add)
  }
  
  # Combine signature results
  all_weights = cbind(weights_snv[,1:ncol(weights_snv)-1],weights_mnv[,1:ncol(weights_mnv)-1], events_id)
  all_weights$sample = id
  return(all_weights)
}

# Bootstrapping function
replicateSignatureChanges <- function(multi_list, time){
  
  # Events list should have all events, even if 0
  all_events = c(rownames(snv_comp), rownames(mnv_comp), "ID1", "ID2", "ID13", "ID8")
  multi_list = lapply(multi_list, function(x){
    add = rep(0, length(all_events[which(!all_events %in% names(x))]))
    names(add) = all_events[which(!all_events %in% names(x))]
    x = c(x, add)
    x[names(x) %in% c("ID1", "ID2", "ID13", "ID8")] = x[names(x) %in% c("ID1", "ID2", "ID13", "ID8")] + 0.001 
    x
  })
  
  a_List_resample = lapply(multi_list, function(x){ 
    rs = t(rmultinom(1, size=sum(x), prob=x)) 
    x = as.numeric(rs)
    names(x) = colnames(rs)
    x
  })
  
  # get weights from new mutations 
  weights_resample = extractSigs(events=a_List_resample)
  
  weights_resample$SBS2.13 = weights_resample$SBS2 + weights_resample$SBS13
  weights_resample$SBS7 = weights_resample$SBS7a + weights_resample$SBS7b + weights_resample$SBS7c + weights_resample$SBS7d
  weights_resample$SBS10 = weights_resample$SBS10a + weights_resample$SBS10b + weights_resample$SBS10c + weights_resample$SBS10d
  weights_resample$SBS17 = weights_resample$SBS17a + weights_resample$SBS17b 
  weights_resample$SBS22 = weights_resample$SBS22a + weights_resample$SBS22b
  weights_resample$SBS40 = weights_resample$SBS40a + weights_resample$SBS40b + weights_resample$SBS40c
  weights_resample$SBS6.14.15.20.21.26.44 = weights_resample$SBS6 + weights_resample$SBS14 + weights_resample$SBS15 + weights_resample$SBS20 + weights_resample$SBS21 + weights_resample$SBS26 + weights_resample$SBS44
  weights_resample = subset(weights_resample[,-which(names(weights_resample) %in% c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS10a", 
                                                               "SBS10b","SBS10c","SBS10d" ,"SBS17a", "SBS17b", "SBS22a","SBS22b","SBS40a","SBS40b","SBS40c","SBS2","SBS13", 
                                                               "SBS6", "SBS14", "SBS15", "SBS20", "SBS26","SBS21", "SBS44"))])
  
  # make proportional, using same totals as for actual changes
  weights_resample$time_total = timed_muts[match(rownames(weights_resample), names(timed_muts))]
  
  clonal_parts <- intersect(rownames(weights_resample), c("clonal [early]", "clonal [late]", "clonal [NA]"))
  if (length(clonal_parts) > 0) {
    clonal_sum <- colSums(weights_resample[clonal_parts, sig_cols, drop = FALSE], na.rm = TRUE)
    weights_resample <- rbind(weights_resample, setNames(rep(NA_real_, ncol(weights_resample)), colnames(weights_resample)))
    rownames(weights_resample)[nrow(weights_resample)] <- "clonal"
    weights_resample["clonal", sig_cols] <- clonal_sum
    
    # time_total / n_unassigned
    weights_resample["clonal", "time_total"] <- sum(weights_resample[clonal_parts, "time_total"], na.rm = TRUE)
    weights_resample["clonal", "n_unassigned"] <- weights_resample["clonal", "time_total"] -
      sum(weights_resample["clonal", sig_cols], na.rm = TRUE)
  }
  
  sig_cols <- grep("^SBS|^DBS|^ID", colnames(weights_resample), value = TRUE)
  weights_resample$n_unassigned = weights_resample$time_total - rowSums(weights_resample[,sig_cols])
  weights_resample[,sig_cols][,which(colSums(weights_resample[,sig_cols]) == 0)] = NA
  prop_cols <- c(sig_cols, "n_unassigned")
  weights_resample[,prop_cols] = weights_resample[,prop_cols]/rowSums(weights_resample[,prop_cols], na.rm=TRUE)
  weights_resample$time = rownames(weights_resample)
  
  weights_resample$n_unassigned[weights_resample$n_unassigned < 0] = 0
  weights_resample[rowSums(weights_resample[,sig_cols], na.rm=TRUE) > 1,prop_cols] = weights_resample[rowSums(weights_resample[,sig_cols], na.rm=TRUE) > 1,prop_cols]/rowSums(weights_resample[rowSums(weights_resample[,sig_cols], na.rm=TRUE) > 1,prop_cols], na.rm=TRUE)
  
  # Calculate changes from early-late and clonal-subclonal
  weights_resample$sample <- id
  if (!"time" %in% colnames(weights_resample)) {
    weights_resample$time <- rownames(weights_resample)
  }
  
  ## Columns to participate in melt: sample, time, all signature columns + n_unassigned
  sig_cols  <- grep("^(SBS|DBS|ID)", colnames(weights_resample), value = TRUE)
  keep_cols <- c("sample", "time", sig_cols, "n_unassigned")
  keep_cols <- intersect(keep_cols, colnames(weights_resample))
  
  ## Long table
  aw_long <- reshape2::melt(
    weights_resample[, keep_cols, drop = FALSE],
    id.vars = c("sample", "time"),
    variable.name = "signature",
    value.name = "weight"
  )
  
  ## Wide table: with signature as row id, expanded into columns by time
  all_weights_long_resample <- reshape(
    aw_long,
    timevar   = "time",
    idvar     = c("sample", "signature"),
    direction = "wide"
  )
  colnames(all_weights_long_resample)[3:ncol(all_weights_long_resample)] = weights_resample$time
  
  # Add pseudocounts and rescale
  all_weights_long_resample[,3:ncol(all_weights_long_resample)] = all_weights_long_resample[,3:ncol(all_weights_long_resample)] + 0.001
  all_weights_long_resample[,3:ncol(all_weights_long_resample)] = sweep(all_weights_long_resample[,3:ncol(all_weights_long_resample)],2,colSums(all_weights_long_resample[,3:ncol(all_weights_long_resample)], na.rm=TRUE),`/`)
  
  
  if (time=="el"){
    # If there are early and late columns, calculate early and late changes
    all_weights_long_resample$early_late = log2((all_weights_long_resample$`clonal [late]`/(1 - all_weights_long_resample$`clonal [late]`)) / (all_weights_long_resample$`clonal [early]`/(1 - all_weights_long_resample$`clonal [early]`)))
    changes = all_weights_long_resample$early_late
    names(changes) = all_weights_long_resample$signature
    changes = as.data.frame.matrix(t(as.data.frame(changes)))
  } else if (time=="cs"){
    # Same for clonal and subclonal
    all_weights_long_resample$clonal_subclonal = log2((all_weights_long_resample$subclonal/(1 - all_weights_long_resample$subclonal)) / (all_weights_long_resample$clonal/(1 - all_weights_long_resample$clonal)))
    changes = all_weights_long_resample$clonal_subclonal
    names(changes) = all_weights_long_resample$signature
    changes = as.data.frame.matrix(t(as.data.frame(changes)))
  }
  
  return(changes)
}




# Read in signatures
snv_comp = read.delim(sbs_path, check.names = FALSE, stringsAsFactors = FALSE)
mnv_comp = read.delim(dbs_path, check.names = FALSE, stringsAsFactors = FALSE)

# Change formatting
## ---- SBS: Generate 96 rownames ----
if ("Type" %in% names(snv_comp) && grepl("\\[", snv_comp$Type[1])) {
  rownames(snv_comp) <- snv_comp$Type
} else if (all(c("Substitution.Type","Trinucleotide") %in% names(snv_comp))) {
  # Another common format (PCAWG style)
  sbs_type <- paste0(substr(snv_comp$Trinucleotide, 1, 1),
                     "[", snv_comp$Substitution.Type, "]",
                     substr(snv_comp$Trinucleotide, 3, 3))
  rownames(snv_comp) <- sbs_type
} else {
  stop("SBS signatures: Unrecognized column naming; please check if the file contains Type or (Substitution.Type, Trinucleotide)")
}

# Remove non-numeric columns, resulting in a matrix of 96 x Nsig.
snv_num_cols <- which(sapply(snv_comp, is.numeric))
snv_comp <- as.matrix(snv_comp[, snv_num_cols, drop = FALSE])


rownames(mnv_comp) = mnv_comp$Type
mnv_comp$Type = NULL
MNV_LEVELS = rownames(mnv_comp)

# Get the multinomials
events = getSNVs_MNVs_IDs(vcf_snv_mnv, vcf_indel, t_df)

snv_samples <- colnames(snv_comp)
# Extract the signatures
n_weights = extractSigs(events=events)

# Combine the signatures here
n_weights$SBS2.13 = n_weights$SBS2 + n_weights$SBS13
n_weights$SBS7 = n_weights$SBS7a + n_weights$SBS7b + n_weights$SBS7c + n_weights$SBS7d
n_weights$SBS10 = n_weights$SBS10a + n_weights$SBS10b + n_weights$SBS10c + n_weights$SBS10d
n_weights$SBS17 = n_weights$SBS17a + n_weights$SBS17b 
n_weights$SBS22 = n_weights$SBS22a + n_weights$SBS22b
n_weights$SBS40 = n_weights$SBS40a + n_weights$SBS40b + n_weights$SBS40c
n_weights$SBS6.14.15.20.21.26.44 = n_weights$SBS6 + n_weights$SBS14 + n_weights$SBS15 + n_weights$SBS20 + n_weights$SBS21 + n_weights$SBS26 + n_weights$SBS44
n_weights = subset(n_weights[,-which(names(n_weights) %in% c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS10a", 
                                                             "SBS10b","SBS10c","SBS10d" ,"SBS17a", "SBS17b", "SBS22a","SBS22b","SBS40a","SBS40b","SBS40c","SBS2","SBS13", 
                                                             "SBS6", "SBS14", "SBS15", "SBS20", "SBS26","SBS21", "SBS44"))])
#n_weights = n_weights[,c(1,65,2:4,69,66,5:6,67,7:9,68,10:64)]

# Add a column of the number of mutations per time frame
timed_muts = t_df
timed_muts = subset(timed_muts, mut_type != "SV")
timed_muts = table(timed_muts$timing)
timed_muts = c(timed_muts, clonal=sum(timed_muts[grep("clonal ", names(timed_muts))]), total=sum(timed_muts))
n_weights$time_total = timed_muts[match(rownames(n_weights), names(timed_muts))]
sig_cols <- grep("^SBS|^DBS|^ID", colnames(n_weights), value = TRUE)
clonal_parts <- intersect(rownames(n_weights), c("clonal [early]", "clonal [late]", "clonal [NA]"))
if (length(clonal_parts) > 0) {
  clonal_sum <- colSums(n_weights[clonal_parts, sig_cols, drop = FALSE], na.rm = TRUE)
  n_weights <- rbind(n_weights, setNames(rep(NA_real_, ncol(n_weights)), colnames(n_weights)))
  rownames(n_weights)[nrow(n_weights)] <- "clonal"
  n_weights["clonal", sig_cols] <- clonal_sum
  
  # time_total / n_unassigned
  n_weights["clonal", "time_total"] <- sum(n_weights[clonal_parts, "time_total"], na.rm = TRUE)
  n_weights["clonal", "n_unassigned"] <- n_weights["clonal", "time_total"] -
    sum(n_weights["clonal", sig_cols], na.rm = TRUE)
}
n_weights$n_unassigned = n_weights$time_total - rowSums(n_weights[, sig_cols])

# Make proportional
all_weights = n_weights

# Where all values are 0, make NA
all_weights[,sig_cols][,which(colSums(all_weights[,sig_cols]) == 0)] = NA
all_weights$time = rownames(all_weights)

# Make proportional
prop_cols <- c(sig_cols, "n_unassigned")
all_weights[,prop_cols] = all_weights[,prop_cols]/all_weights$time_total

# If more mutations are assigned than are actually in the sample, rescale so sum to 1
all_weights$n_unassigned[all_weights$n_unassigned < 0] = 0
all_weights[rowSums(all_weights[,sig_cols], na.rm=TRUE) > 1, prop_cols] = all_weights[rowSums(all_weights[,sig_cols], na.rm=TRUE) > 1, prop_cols]/rowSums(all_weights[rowSums(all_weights[,sig_cols], na.rm=TRUE) > 1, prop_cols], na.rm=TRUE)


# Calculate changes from early-late and clonal-subclonal
all_weights$sample <- id
if (!"time" %in% colnames(all_weights)) {
  all_weights$time <- rownames(all_weights)
}
## Columns to participate in melt: sample, time, all signature columns + n_unassigned
sig_cols  <- grep("^(SBS|DBS|ID)", colnames(all_weights), value = TRUE)
keep_cols <- c("sample", "time", sig_cols, "n_unassigned")
keep_cols <- intersect(keep_cols, colnames(all_weights))  # 以防万一

## long table
aw_long <- reshape2::melt(
  all_weights[, keep_cols, drop = FALSE],
  id.vars = c("sample", "time"),
  variable.name = "signature",
  value.name = "weight"
)

## wide table: with signature as row id, expanded into columns by time
all_weights_long <- reshape(
  aw_long,
  timevar   = "time",
  idvar     = c("sample", "signature"),
  direction = "wide"
)
colnames(all_weights_long)[3:ncol(all_weights_long)] = all_weights$time

# Add pseudocounts and rescale
all_weights_long[,3:ncol(all_weights_long)] = all_weights_long[,3:ncol(all_weights_long)] + 0.001
all_weights_long[,3:ncol(all_weights_long)] = sweep(all_weights_long[,3:ncol(all_weights_long)],2,colSums(all_weights_long[,3:ncol(all_weights_long)], na.rm=TRUE),`/`)

# If there are early and late columns, calculate early and late changes
if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
  all_weights_long$early_late = log2((all_weights_long$`clonal [late]`/(1 - all_weights_long$`clonal [late]`)) / (all_weights_long$`clonal [early]`/(1 - all_weights_long$`clonal [early]`)))
}

# Same for clonal and subclonal
if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
  all_weights_long$clonal_subclonal = log2((all_weights_long$subclonal/(1 - all_weights_long$subclonal)) / (all_weights_long$clonal/(1 - all_weights_long$clonal)))
}


# Bootstrapping for early-late and clonal-subclonal

if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
  changes_1000_el = as.data.frame.matrix(t(replicate(500, replicateSignatureChanges(multi_list=events, time="el"), simplify="vector")))
}

if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
  changes_1000_cs = as.data.frame.matrix(t(replicate(500, replicateSignatureChanges(multi_list=events, time="cs"), simplify="vector")))
}
  
# Get confidence intervals for early-late signature changes
all_weights_long = all_weights_long[!is.na(rowSums(all_weights_long[,3:ncol(all_weights_long)])),]

CIs = apply(all_weights_long, 1, function(x) {
  
  sig = x[2] # which signature?
  
  # Early/late confidence intervals
  if ("clonal [early]" %in% colnames(all_weights_long) & "clonal [late]" %in% colnames(all_weights_long)){
    el_reps = unlist(changes_1000_el[,sig])
    el_lower = quantile(as.numeric(el_reps), 0.025, na.rm=TRUE)
    el_upper = quantile(as.numeric(el_reps), 0.975, na.rm=TRUE)
  } else { 
    el_lower=NA
    el_upper=NA
  }
 

  # Clonal/subclonal
  if ("clonal" %in% colnames(all_weights_long) & "subclonal" %in% colnames(all_weights_long)){
    cs_reps = unlist(changes_1000_cs[,sig])
    cs_lower = quantile(as.numeric(cs_reps), 0.025, na.rm=TRUE)
    cs_upper = quantile(as.numeric(cs_reps), 0.975, na.rm=TRUE)  
  } else {
    cs_lower=NA
    cs_upper=NA
  }
  
  c(el_lower, el_upper, cs_lower, cs_upper)
  
})



all_weights_long$el_lCI = CIs[1,]
all_weights_long$el_uCI = CIs[2,]
all_weights_long$cs_lCI = CIs[3,]
all_weights_long$cs_uCI = CIs[4,]

# Add n to all_weights, remove if both 0
n_weights$time = rownames(n_weights)
n_weights$sample <- id

n_weights_m <- reshape2::melt(
  n_weights[, keep_cols, drop = FALSE],
  id.vars = c("sample", "time"),
  variable.name = "signature",
  value.name = "weight"
)

n_weights_m = reshape(n_weights_m,
                      timevar="time",
                      idvar=c("sample", "signature"),
                      direction="wide")
n_prop_weights = merge(all_weights_long, n_weights_m)

# Save all results 
write.table(all_weights_long, paste0(id, "_sig_weights.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Combine bootstrapping results
if (exists("changes_1000_cs")){
  changes_1000_cs = apply(changes_1000_cs, 2, unlist)
  mode(changes_1000_cs) = "numeric"
  changes_1000_cs = data.frame(changes_1000_cs)
  changes_1000_cs$time = "clonal_subclonal"
} 

if (exists("changes_1000_el")){
  changes_1000_el = apply(changes_1000_el, 2, unlist)
  mode(changes_1000_el) = "numeric"
  changes_1000_el = data.frame(changes_1000_el)
  changes_1000_el$time = "early_late"
}

#if (exists("changes_1000_el") & exists("changes_1000_cs")){
  
#  all_bootstraps = rbind(changes_1000_cs, changes_1000_el)
#  write.table(all_bootstraps, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", row.names=TRUE, quote=FALSE)

#} else if(exists("changes_1000_cs") & !exists("changes_1000_el")){
#    write.table(changes_1000_cs, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", row.names=TRUE, quote=FALSE)
  
#} else if(exists("changes_1000_el") & !exists("changes_1000_cs")){
#    write.table(changes_1000_el, file=gzfile(paste0(id, "_sig_change_bootstraps.txt.gz")), sep="\t", row.names=TRUE, quote=FALSE)
  
#}

###########################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# ---- Define the standard order of COSMIC 96 channels ----
make_cosmic96_levels <- function() {
  subs  <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  bases <- c("A","C","G","T")
  as.vector(unlist(lapply(subs, function(s)
    outer(bases, bases, FUN = function(b5,b3) sprintf("%s[%s]%s", b5, s, b3)))))
}

cosmic96_levels <- make_cosmic96_levels()
as96 <- rownames(snv_comp)
# Here is an alignment made
as96 <- intersect(cosmic96_levels, as96)

# ------- Translate events into a long table (count) and label 6 types of replacements -------
make_counts_df <- function(events, as96, timings = c("clonal [early]", "clonal [late]")) {
  timings <- intersect(timings, names(events))
  stopifnot(length(timings) >= 1)
  
  mat <- sapply(events[timings], function(v){
    x <- integer(length(as96)); names(x) <- as96
    v <- v[intersect(names(v), as96)]
    x[names(v)] <- as.integer(v)
    x
  })
  if (is.null(dim(mat))) mat <- matrix(mat, ncol=1, dimnames=list(as96, timings))
  df <- melt(mat, varnames = c("Type96","Timing"), value.name = "Count")
  
  df$Substitution <- sub("^.*\\[([ACGT]>[ACGT])\\].*$", "\\1", df$Type96)
  df$Type96 <- factor(df$Type96, levels = as96)
  df$Substitution <- factor(df$Substitution, levels = c("C>A","C>G","C>T","T>A","T>C","T>G"))
  df
}

sub6_cols <- c(
  "C>A"="#3db7e9",  # blue
  "C>G"="#000000",  # black
  "C>T"="#e84a5f",  # red
  "T>A"="#bfbfbf",  # grey
  "T>C"="#65c36f",  # green
  "T>G"="#f4a6c1"   # pink
)

# ------- Plot a side-by-side 96-channel digital spectrum (y-axis free) -------
plot_96_pair_counts <- function(events, as96,
                                timings = c("clonal [early]","clonal [late]"),
                                title = NULL) {
  df <- make_counts_df(events, as96, timings)
  ggplot(df, aes(x = Type96, y = Count, fill = Substitution)) +
    geom_col(width = 0.9) +
    facet_grid(~ Timing, scales = "free_y") +
    scale_fill_manual(values = sub6_cols, drop = FALSE) +
    labs(x = "96-channel context", y = "SNVs", title = title) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      strip.background = element_rect(fill = "white", colour = "grey80")
    )
}

# -------Annotation of signature proportion at each time layer -------
annotate_sig_perc <- function(p, prop_df, timings = c("clonal [early]","clonal [late]"),
                              sig_names = c("SBS7","SBS38")) {
  if (missing(prop_df) || is.null(prop_df)) return(p)
  ann <- lapply(intersect(timings, unique(prop_df$time)), function(ti){
    row <- prop_df[prop_df$time == ti, , drop = FALSE][1, , drop = FALSE]
    av <- row[ , intersect(sig_names, colnames(prop_df)), drop = FALSE]
    if (ncol(av) == 0) return(NULL)
    txt <- paste0(
      paste0(round(100*as.numeric(av[1,]), 0), "% ", colnames(av)),
      collapse = "\n"
    )
    data.frame(Timing = ti, label = txt, stringsAsFactors = FALSE)
  })
  ann <- do.call(rbind, ann)
  if (is.null(ann)) return(p)
  
  p + geom_text(
    data = ann, aes(x = 5, y = Inf, label = label),  # x=5 Place text at the top left; adjustable
    inherit.aes = FALSE, vjust = 1.05, hjust = 0,
    size = 3.3
  )
}

# ======== Early vs Late ========
p_el <- plot_96_pair_counts(
  events, as96,
  timings = c("clonal [early]", "clonal [late]"),
  title = paste0(id, " — Early clonal vs Late clonal (96-channel counts)")
)

# ======== Clonal vs Subclonal ========

sum_counts <- function(...) {
  xs <- list(...)
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (length(xs) == 0) return(setNames(numeric(0), character(0)))
  
  all_names <- Reduce(union, lapply(xs, names))
  aligned <- lapply(xs, function(v) {
    out <- setNames(numeric(length(all_names)), all_names)
    out[names(v)] <- as.numeric(v)
    out
  })
  Reduce(`+`, aligned)
}

# Layers to be merged
parts <- intersect(c("clonal [early]", "clonal [late]", "clonal [NA]"), names(events))

# If these layers exist, generate new clonal layers.
if (length(parts) > 0) {
  events[["clonal"]] <- do.call(sum_counts, events[parts])
}

p_cs <- plot_96_pair_counts(
  events, as96,
  timings = c("clonal", "subclonal"),
  title = paste0(id, " — Clonal vs Subclonal (96-channel counts)")
)



## ========= Help function: Extract the maximum positive/negative change, and obtain the percentage before/after. =========
get_top_changes_with_props <- function(df, change_col, before_col, after_col) {
  stopifnot(all(c(change_col, before_col, after_col) %in% colnames(df)))
  d <- df[, c("signature", change_col, before_col, after_col)]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  pos <- d[which.max(d[[change_col]]), , drop = FALSE]
  neg <- d[which.min(d[[change_col]]), , drop = FALSE]
  list(
    pos = list(sig = pos$signature,
               change = as.numeric(pos[[change_col]]),
               before = as.numeric(pos[[before_col]]),
               after  = as.numeric(pos[[after_col]])),
    neg = list(sig = neg$signature,
               change = as.numeric(neg[[change_col]]),
               before = as.numeric(neg[[before_col]]),
               after  = as.numeric(neg[[after_col]]))
  )
}

## ========= Data required for generating geom_text (including percentage from front to back and ×FC) =========
make_ann_with_props <- function(top_list, panel_name, x = 6, y = Inf) {
  if (is.null(top_list)) return(NULL)
  fmt_pct <- function(x) paste0(sprintf("%0.1f", 100 * x), "%")
  fmt_line <- function(tag, x) {
    paste0(
      tag, " ", x$sig, ": ",
      fmt_pct(x$before), " → ", fmt_pct(x$after),
      "  (", sprintf("%+0.2f", x$change), ", ×", sprintf("%0.2f", 2^x$change), ")"
    )
  }
  data.frame(
    Timing = panel_name,
    label  = paste0(
      "↑ ", fmt_line("", top_list$pos), "\n",
      "↓ ", fmt_line("", top_list$neg)
    ),
    x = x, y = y,
    stringsAsFactors = FALSE
  )
}

## ========= Early vs Late  =========
ann_el <- NULL
if (all(c("clonal [early]","clonal [late]") %in% colnames(all_weights_long)) &&
    "early_late" %in% colnames(all_weights_long)) {
  tl_el <- get_top_changes_with_props(
    all_weights_long,
    change_col = "early_late",
    before_col = "clonal [early]",
    after_col  = "clonal [late]"
  )
  ann_el <- make_ann_with_props(tl_el, panel_name = "clonal [late]", x = 6)
  if (!is.null(ann_el)) {
    p_el <- p_el +
      ggplot2::geom_text(
        data = ann_el,
        mapping = aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        vjust = 1.05, hjust = 0, size = 3.3
      )
  }
}

## ========= Clonal vs Subclonal =========
ann_cs <- NULL
if (all(c("clonal","subclonal") %in% colnames(all_weights_long)) &&
    "clonal_subclonal" %in% colnames(all_weights_long)) {
  tl_cs <- get_top_changes_with_props(
    all_weights_long,
    change_col = "clonal_subclonal",
    before_col = "clonal",
    after_col  = "subclonal"
  )
  ann_cs <- make_ann_with_props(tl_cs, panel_name = "subclonal", x = 6)
  if (!is.null(ann_cs)) {
    p_cs <- p_cs +
      ggplot2::geom_text(
        data = ann_cs,
        mapping = aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        vjust = 1.05, hjust = 0, size = 3.3
      )
  }
}

## save plots
ggsave(paste0(id, "_early_vs_late_96counts.pdf"), p_el, width = 14, height = 4.5)
ggsave(paste0(id, "_clonal_vs_subclonal_96counts.pdf"), p_cs, width = 14, height = 4.5)