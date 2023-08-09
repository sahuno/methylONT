#samuel ahuno
#generate promoter regions for genes; 1kb upstream and 1kb downstream
## code to prepare `gene_promoters_encode1kb` dataset goes here
#library(methylONT)
library(data.table)
library(dplyr)
library(annotables)
library(BSgenome.Hsapiens.UCSC.hg38)
#load("/juno/work/greenbaum/users/ahunos/apps/methylONT/data/gene_promoters_encode1kb.rda")
#gene_promoters_encode1kb

#specify paths to data
path_2_hg38TSS <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/ENCFF493CCB.bed.gz"


dt <- fread(path_2_mm10TSS)
#head(dt)
names(dt) <- c("seqnames", "start", "end", "ensgene.version", "score", "strand")
#dim(dt)
dt[, ensgene:=(gsub("\\..*","",ensgene.version))]
dt[, .N, by = ensgene.version][N>1,] #sanity check, is any transcript annotated more than once?
dt[, .N, by = ensgene][N>1,] #sanity check, is any transcript annotated more than once? yes, there are about 18 transcripts with multiple promoter annotations



# dt[grep("ENSG00000133703", V4),]
# grch38_tx2gene %>% dplyr::filter(str_detect(enstxp, "ENSG00000133703"))
# grch38_tx2gene_dt <- grch38_tx2gene %>% as.data.table()
# dt[grch38_tx2gene_dt, on = "ensgene", nomatch = NULL]

grch38_dt <- grch38 %>% as.data.table()
grch38_dt <- grch38_dt[, .(ensgene, entrez, symbol, biotype, description)]
dt <- dt[grch38_dt, on = "ensgene", nomatch = NULL]
#dt[symbol == "",] # remove genes without names 
#na.omit(dt, "symbol")
dt <- dt[, .(seqnames, start, end, strand, ensgene.version, symbol, biotype, key = paste0(symbol,"_" ,ensgene.version, strand))] #save only needed parts 
setkey(dt, key)

dt <- unique(dt, by = "key") #remove duplicate promoters
# dt[, fullDup := .N > 1, by = key(dt)] #check for full duplicates in the data
# dt_keyed =  dt_keyed[fullDup == ,]

#add to utils scripts
message("adding transcription start site")
dt[,`:=`(transcription_start_site = fcase(strand == "+", start, strand == "-", end))] #set transcription start site

dt[,.N, by = biotype] #check biotype

gr_tss_encode <- GenomicRanges::makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
gene_promoters_encode1kb <- GenomicRanges::promoters(gr_tss_encode, upstream = 1000, downstream = 1000)

#get num_GCs in gene promoters 
gene_promoters_encode1kb_Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, gene_promoters_encode1kb)
gene_promoters_encode1kb_numGCs <- Biostrings::vcountPattern("GC", gene_promoters_encode1kb_Seq)
gene_promoters_encode1kb_numCpGs <- Biostrings::vcountPattern("CG", gene_promoters_encode1kb_Seq)

mcols(gene_promoters_encode1kb)$numGCs <- gene_promoters_encode1kb_numGCs
mcols(gene_promoters_encode1kb)$numCGs <- gene_promoters_encode1kb_numCpGs

gene_promoters_encode1kb_proteinCoding <- gene_promoters_encode1kb[gene_promoters_encode1kb$biotype == "protein_coding"]
#protein_coding
#uncomment to make promoters for encode genes
#usethis::use_data_raw("gene_promoters_encode1kb")
usethis::use_data(gene_promoters_encode1kb, overwrite = TRUE)
usethis::use_data(gene_promoters_encode1kb_proteinCoding, overwrite = TRUE)



##############################################################################
##############################################################################
#create data sets 
ff <- "/work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-044_N/per_read_modified_base_calls.txt"
dt1 <- data.table::fread(ff, nrows = 100)
cmd = glue::glue("tail -n 100 {ff}")
dt2 <- data.table::fread(cmd = cmd, header=TRUE)
per_mod_base_DT <- data.table::rbindlist(list(dt1, dt2))
usethis::use_data(per_mod_base_DT, overwrite = TRUE)

##testing how the reads are organized
# load("/juno/work/greenbaum/users/ahunos/apps/methylONT/data/per_mod_base_DT.rda")
dt <- copy(per_mod_base_DT)
dt2 <- dt[!grep("_",chrm),][, `:=`(exp_mod_log_prob = exp(mod_log_prob), exp_can_log_prob = exp(can_log_prob))][,`:=`(sum_Exp_mod_can_prob=exp_mod_log_prob+exp_can_log_prob, sum_mod_can_prob=mod_log_prob+can_log_prob)]
summary(dt2$sum_mod_can_prob)
summary(dt2$sum_Exp_mod_can_prob)                                    
split(dt2, dt2$mod_base)
###################################################################################




###################################################################################
#############Mouse data promoters and enhancers ####################
###################################################################################
path_2_mm10TSS <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/encode/ATAC-seq/mm10/ENCFF498BEJ.bed.gz"

dt <- fread(path_2_mm10TSS)
names(dt) <- c("seqnames", "start", "end", "ensgene.version", "score", "strand")
dt[, ensgene:=(gsub("\\..*","",ensgene.version))]
#dt[, .N, by = ensgene.version][N>1,] #sanity check, is any transcript annotated more than once?

grcm38_dt <- grcm38 %>% as.data.table()
grcm38_dt <- grcm38_dt[, .(ensgene, entrez, symbol, biotype, description)]
dt <- dt[grcm38_dt, on = "ensgene", nomatch = NULL]

dt <- dt[, .(seqnames, start, end, strand, ensgene.version, symbol, biotype, key = paste0(symbol,"_" ,ensgene.version, strand))] #save only needed parts 
setkey(dt, key)

dt <- unique(dt, by = "key") #remove duplicate promoters
# dt[, fullDup := .N > 1, by = key(dt)] #check for full duplicates in the data
# dt_keyed =  dt_keyed[fullDup == ,]

#add to utils scripts
message("adding transcription start site")
dt[,`:=`(transcription_start_site = fcase(strand == "+", start, strand == "-", end))] #set transcription start site

dt[,.N, by = biotype] #check biotype

gr_tss_encode_mm10 <- GenomicRanges::makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
gene_promoters_encode1kb_mm10 <- GenomicRanges::promoters(gr_tss_encode_mm10, upstream = 1000, downstream = 1000)

#get num_GCs in gene promoters 
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
#head(seqlengths(genome), 23)
gene_promoters_encode1kb_mm10_Seq <- getSeq(genome, gene_promoters_encode1kb_mm10)
gene_promoters_encode1kb_mm10_numGCs <- Biostrings::vcountPattern("GC", gene_promoters_encode1kb_mm10_Seq)
gene_promoters_encode1kb_mm10_numCpGs <- Biostrings::vcountPattern("CG", gene_promoters_encode1kb_mm10_Seq)

mcols(gene_promoters_encode1kb_mm10)$numGCs <- gene_promoters_encode1kb_mm10_numGCs
mcols(gene_promoters_encode1kb_mm10)$numCGs <- gene_promoters_encode1kb_mm10_numCpGs

gene_promoters_encode1kb_mm10_proteinCoding <- gene_promoters_encode1kb_mm10[gene_promoters_encode1kb_mm10$biotype == "protein_coding"]

saveRDS(gene_promoters_encode1kb_mm10_proteinCoding, file = "inst/extdata/gene_promoters_encode1kb_mm10_proteinCoding.rds")
saveRDS(gene_promoters_encode1kb_mm10, file = "inst/extdata/gene_promoters_encode1kb_mm10.rds")

#protein_coding
#uncomment to make promoters for encode genes
#usethis::use_data_raw("gene_promoters_encode1kb")
usethis::use_data(gene_promoters_encode1kb_mm10, overwrite = TRUE)
usethis::use_data(gene_promoters_encode1kb_mm10_proteinCoding, overwrite = TRUE)
