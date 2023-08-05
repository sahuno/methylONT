#' function to export megalodon aggregate mod bases to bismark like coverage file
#' chromosome start end methylation_percentage count_methylated count_unmethylated
#' @param dt data.table object
#' @param file_stem stem of file name to write to
#' @return aggregate statistics of 5mc and 5hmc from `aggregate per read modified bases from megalodon`; 
#' 
#' @export
#' @import data.table
#' @import utils
make_bismark_coverage <- function(dt, file_stem=NULL){
bismark_coverage_dt <- dt[,.(seqnames, end, end, strand, methylation_percentage_5mC_5hmC, count_methylated, count_unmethylated)]  
data.table::fwrite(bismark_coverage_dt, file = paste0(file_stem,"ont_aggregate_5hmC_5mC_bismark.cov.gz"), sep = '\t', col.names = FALSE)
}
utils::globalVariables(c("seqnames", "end", "methylation_percentage_5mC_5hmC", "count_methylated", "strand", "count_unmethylated"))

## @examples 
## make_bismark_coverage(dt=DT, file_stem="dir/sampleXYZ_")