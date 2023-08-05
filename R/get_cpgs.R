#' function to make aggregate statistics from per read modified probabilities from megalodon
#'
#' @param species which species to use, options "human" or "mouse"; implement for mouse
#' @param genome which genome to use, options "hg19" or "hg38"
#' @param offset 0/1-based granges object
#' @return granges object of cpg loci
#' 
#' @export
#' @import AnnotationHub
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import utils
#' @examples 
#' cpg_sites(species="human",genome="hg38", offset=0)
cpg_sites <- function(species="human",genome="hg38", offset=0){
############################
###get cpg loci
############################
  if(species=="human" & genome=="hg38"){
    chr_keep <- c(paste0("chr",c(seq(22), "X","Y"))) #specify chrms to extract
    length_chrsHg38 <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[chr_keep] #get lenghts of chroms of intrest
    cgs <- lapply(chr_keep, function(x) start(Biostrings::matchPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]))) #cpgs loci per chrm
        if(offset==0){
      cpgr_0based <- do.call(c, lapply(seq_along(chr_keep), function(x) GenomicRanges::GRanges(names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[x], IRanges::IRanges(start = cgs[[x]]-1, width = 2)))) 
      GenomeInfoDb::seqlengths(cpgr_0based) <- length_chrsHg38 #set chrom lengths from information from hg19
      return(cpgr_0based)
        }else{
      cpgr_1based <- do.call(c, lapply(seq_along(chr_keep), function(x) GenomicRanges::GRanges(names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[x], IRanges::IRanges(start = cgs[[x]], width = 1))))  #this is 1-based, verified with ensemble browser whhicj is 1-based https://grch37.ensembl.org/Homo_sapiens/Location/View?r=1%3A10469-10469
      GenomeInfoDb::seqlengths(cpgr_1based) <- length_chrsHg38 #set chrom lengths from information from hg19
      return(cpgr_1based)
        }
    
  }
  #else if(species=="human" & genome=="hg19"){
  #   require("BSgenome.Hsapiens.UCSC.hg19")
  #   length_chrsHg19 <- seqlengths(Hsapiens)[chr_keep] #get lenghts of chroms of intrest
  #       if(offset==0){
  #         cpgr_0based <- do.call(c, lapply(seq_along(chr_keep), function(x) GRanges(names(Hsapiens)[x], IRanges(start = cgs[[x]]-1, width = 2)))) 
  #         seqlengths(cpgr_0based) <- length_chrsHg19 #set chrom lengths from information from hg19
  #       }else{
  #         cpgr_1based <- do.call(c, lapply(seq_along(chr_keep), function(x) GRanges(names(Hsapiens)[x], IRanges(start = cgs[[x]], width = 1))))  #this is 1-based, verified with ensemble browser whhicj is 1-based https://grch37.ensembl.org/Homo_sapiens/Location/View?r=1%3A10469-10469
  #         seqlengths(cpgr_1based) <- length_chrsHg19 #set chrom lengths from information from hg19
  #   }
  # }
}

#test_cpg_default <- cpg_sites()

#    cgs_counts <- lapply(chr_keep, function(x) start(Biostrings::matchPattern("CG", BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]))) #cpgs loci per chrm

#countPattern(dnaseq, Scerevisiae$chrI)
# BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]]