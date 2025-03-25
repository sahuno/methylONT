#library(data.table)
#read multiple mpileup bed files
multi_pileup <- function(files_paths){
data_list <- lapply(paths_bed, data.table::fread, col.names = c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "Freq_5mCpG", "mod", "canon", "other", "del", "fail", "diff", "nocall")) 

}


