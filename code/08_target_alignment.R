# Load required libraries
library(biomaRt)
library(Biostrings)
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# Load validated interactions
interactions <- fread("validated_interactions.csv")
setnames(interactions, c("miRNA","Target_Gene"))

# Fetch canonical 3'UTRs
mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "https://www.ensembl.org")
attrs <- c("hgnc_symbol","ensembl_transcript_id","chromosome_name","strand",
           "transcript_start","transcript_end","3utr","transcript_is_canonical")
target_genes <- sort(unique(interactions$Target_Gene))
batch_size <- 300
chunks <- split(target_genes, ceiling(seq_along(target_genes)/batch_size))
utr_list <- list()
for (i in seq_along(chunks)) {
  utr_list[[i]] <- getBM(attributes = attrs,
                         filters = "hgnc_symbol",
                         values = chunks[[i]], mart = mart)
}
utr_data <- as.data.table(rbindlist(utr_list, use.names=TRUE, fill=TRUE))
utr_data <- utr_data[transcript_is_canonical == 1]
setnames(utr_data, names(utr_data), sub("^X3utr$","3utr",names(utr_data)))
utr_data <- utr_data[!is.na(`3utr`) & `3utr` != "" & `3utr` != "Sequence unavailable"]
utr_data[, utr_len := nchar(`3utr`)]
utr_data <- utr_data[chromosome_name %in% c(as.character(1:22),"X","Y","MT","M")]
setorder(utr_data,hgnc_symbol,-utr_len)
utr_data <- utr_data[!duplicated(hgnc_symbol)]
utr_data[, strand := as.integer(strand)]
utr_data[, `:=`(transcript_start=as.integer(transcript_start),
                transcript_end=as.integer(transcript_end))]

# Build UTR DNAStringSet
utr_set <- DNAStringSet(utr_data$`3utr`)
names(utr_set) <- utr_data$hgnc_symbol

# Load miRNAs
mirna_raw <- readRNAStringSet("/mature_mirnas.fa", format="fasta")
mirna_dna <- DNAStringSet(RNAStringSet(mirna_raw))
names(mirna_dna) <- sub(" .*","",trimws(names(mirna_dna)))
mirna_dna <- mirna_dna[grepl("^hsa-",names(mirna_dna))]
mirna_dna <- mirna_dna[names(mirna_dna) %in% interactions$miRNA]

# Build seed dictionaries
seed6   <- DNAStringSet(lapply(mirna_dna,function(s) subseq(s,2,7)))
seed7m8 <- DNAStringSet(lapply(mirna_dna,function(s) subseq(s,2,8)))
seed7A1 <- DNAStringSet(lapply(mirna_dna,function(s) subseq(s,1,7)))
seed8   <- DNAStringSet(lapply(mirna_dna,function(s) subseq(s,1,8)))
dict6   <- PDict(reverseComplement(seed6))
dict7m8 <- PDict(reverseComplement(seed7m8))
dict7A1 <- PDict(reverseComplement(seed7A1))
dict8   <- PDict(reverseComplement(seed8))

# Scan UTRs
res_list <- list()
for (i in seq_along(utr_set)) {
  gene <- names(utr_set)[i]
  subj <- utr_set[[i]]
  hits <- data.table()
  h6 <- matchPDict(dict6,subj)
  if (any(lengths(h6)>0)) hits <- rbind(hits,data.table(Target_Gene=gene,
                                                        miRNA=rep(names(seed6),lengths(h6)),
                                                        bs_start=start(unlist(h6)),
                                                        bs_end=end(unlist(h6)),
                                                        seed_class="6mer"))
  h7m8 <- matchPDict(dict7m8,subj)
  if (any(lengths(h7m8)>0)) hits <- rbind(hits,data.table(Target_Gene=gene,
                                                          miRNA=rep(names(seed7m8),lengths(h7m8)),
                                                          bs_start=start(unlist(h7m8)),
                                                          bs_end=end(unlist(h7m8)),
                                                          seed_class="7mer-m8"))
  h7A1 <- matchPDict(dict7A1,subj)
  if (any(lengths(h7A1)>0)) hits <- rbind(hits,data.table(Target_Gene=gene,
                                                          miRNA=rep(names(seed7A1),lengths(h7A1)),
                                                          bs_start=start(unlist(h7A1)),
                                                          bs_end=end(unlist(h7A1)),
                                                          seed_class="7mer-A1"))
  h8 <- matchPDict(dict8,subj)
  if (any(lengths(h8)>0)) hits <- rbind(hits,data.table(Target_Gene=gene,
                                                        miRNA=rep(names(seed8),lengths(h8)),
                                                        bs_start=start(unlist(h8)),
                                                        bs_end=end(unlist(h8)),
                                                        seed_class="8mer"))
  if (nrow(hits)>0) res_list[[i]] <- hits
}
all_hits <- rbindlist(res_list,use.names=TRUE,fill=TRUE)

# Keep only validated pairs
all_hits[, gene_miRNA := paste(Target_Gene,miRNA,sep="_")]
interactions[, gene_miRNA := paste(Target_Gene,miRNA,sep="_")]
hits <- all_hits[gene_miRNA %in% interactions$gene_miRNA]
hits[, gene_miRNA := NULL]

# Merge transcript info
hits <- merge(hits,
              utr_data[,.(Target_Gene=hgnc_symbol,
                          ensembl_transcript_id,
                          transcript_length,
                          utr_len,
                          `3utr_start_tx`=transcript_start)],
              by="Target_Gene",all.x=TRUE)

# Map to genomic coordinates
utr_by_tx <- threeUTRsByTranscript(EnsDb.Hsapiens.v86)
strip_ver <- function(x) sub("\\..*$","",x)
names(utr_by_tx) <- strip_ver(names(utr_by_tx))
hits[, ensembl_transcript_id := strip_ver(ensembl_transcript_id)]
hits <- hits[!is.na(ensembl_transcript_id) & nzchar(ensembl_transcript_id)]
hits_split <- split(hits,by="ensembl_transcript_id",keep.by=TRUE)

mapped_list <- list()
for (i in seq_along(hits_split)) {
  tx_hits <- hits_split[[i]]
  tx_id <- unique(tx_hits$ensembl_transcript_id)[1]
  if (!(tx_id %in% names(utr_by_tx))) next
  utr_len_tx <- sum(width(utr_by_tx[[tx_id]]))
  tx_hits <- tx_hits[bs_start>=1L & bs_end<=utr_len_tx]
  if (!nrow(tx_hits)) next
  gr_tx <- GRanges(seqnames=Rle(rep(tx_id,nrow(tx_hits))),
                   ranges=IRanges(start=tx_hits$bs_start,end=tx_hits$bs_end))
  mapped <- mapFromTranscripts(gr_tx,utr_by_tx)
  if (length(mapped)==0L) next
  gr_dt <- as.data.table(mapped)
  xh <- mcols(mapped)$xHits
  thits <- cbind(tx_hits[xh],
                 chromosome=as.character(gr_dt$seqnames),
                 genomic_start=gr_dt$start,
                 genomic_end=gr_dt$end,
                 target_strand=as.character(gr_dt$strand))
  mapped_list[[i]] <- thits
}
hits_expanded <- unique(rbindlist(mapped_list,use.names=TRUE,fill=TRUE))

# Export BED and CSV
hits_expanded[, chromosome := ifelse(grepl("^chr",chromosome),chromosome,paste0("chr",chromosome))]
hits_expanded[, target_start := genomic_start-1L]
hits_expanded[, target_end := genomic_end]
hits_expanded[, name := paste(miRNA,Target_Gene,chromosome,target_start,target_end,target_strand,sep="_")]

out <- hits_expanded[,.(chromosome,target_start,target_end,
                        target_strand,name,seed_class,
                        ensembl_transcript_id)]

fwrite(out,"mirna_pairs_genomic.bed",sep="\t",col.names=TRUE,quote=FALSE)
