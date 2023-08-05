#' gene promoters from encode databasae, defined as 1000 bp of transcription start site
#' @format ## `gene_promoters_encode1kb`
#' A GenomicRanges object
#' \describe{
#'  \item{seqnames}{chromosome name}
#'   \item{ranges}{genomic start and end position}
#'  \item{strand}{strand of the gene}
#' \item{ensgene.version}{ensembl gene id with version}
#' \item{ensgene}{ensembl gene id}
#' \item{entrez}{entrez id}
#' \item{biotype}{gene biotype}
#' \item{description}{gene description}
#' \item{symbol}{gene name}
#'}
#' @source <https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz>
"gene_promoters_encode1kb"
