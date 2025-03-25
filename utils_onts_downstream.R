
# hmm, tired of writing functions all the time
library(data.table)
library("fst")

# args(read_fst)
##################################
#read multiple pileups tables with option t
multi_mpileup <- function(paths, colnames){
    # default is mpileup from modkit
    if(null(colnames)){
    ls_data <- lapply(paths, fread, col.names = c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "Freq_5mCpG", "mod", "canon", "other", "del", "fail", "diff", "nocall")) 
    }else{
        ls_data <- lapply(paths_bed, fread, col.names = colnames) 
    }
}

multi_mpileup_fst <- function(paths, colnames){
    # default is mpileup from modkit
    if(null(colnames)){
    ls_data <- lapply(paths, function(x){
        fst::read_fst(x, columns = c("chrom", "start", "end", "name", "score", "strand", "tstart", "tend", "color", "coverage", "Freq_5mCpG", "mod", "canon", "other", "del", "fail", "diff", "nocall")) 
    } ) 
    }else{
        ls_data <- lapply(paths, function(x){
        fst::read_fst(x, columns = colnames) 
    }) 
    }
}

################
#sample name is after first `_`
getSampleNames <- function(x){
    print("pass")
}

########################
##plotting functions


## get sampels in snakemkae pipeline