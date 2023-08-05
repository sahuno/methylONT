#' function to make aggregate statistics from per read modified probabilities from megalodon
#'
#' @param paths_aggregate_mod_bases vector of paths to megalodon aggregate mod bases
#' @param methyl_percentage_threshold minimum percentage_methylation either 5mc or 5hmC to be called methylated otherwise unmethylated
#' @param methyl_numb_reads_threshold minimum number of reads to be called methylated otherwise unmethylated
#' @return data.table object representing coverage of aggregate 5mc and 5hmc from megalodon, methylation percentage & state (5hmc, 5mc, U)
#' @export
#' @import data.table
#' @import utils
#' @examples 
#' #read_megalodon_aggregate_mod_bases(paths_aggregate_mod_bases=c("dir/5hmc.bed", "dir/5mc.bed"), methyl_percentage_threshold=0.5, methyl_numb_reads_threshold=5)

# .datatable.aware <- TRUE
read_megalodon_aggregate_mod_bases <- function(paths_aggregate_mod_bases, methyl_percentage_threshold, methyl_numb_reads_threshold){
 #print(paths_aggregate_mod_bases) 
  #we assume this is a list of 2 file 
  read_megalodon_aggregate_mod_bases_from_file <- function(x){
    #require(data.table)
    dt <- data.table::fread(x, header = FALSE, drop = c(4,7,8,9,10))
    names(dt) <- c("seqnames", "start", "end", "number_reads","strand","methylation_percentage")
    return(dt)
    #droppinig columns: c(4,7,8,9,10) from megalodon data
  }
  
  message(paste0("reading ", paths_aggregate_mod_bases, "...\n"))
  ls_5mC_5hmC <- lapply(paths_aggregate_mod_bases, read_megalodon_aggregate_mod_bases_from_file)
  names(ls_5mC_5hmC) <- basename(paths_aggregate_mod_bases) #rename headers
  ls_5mC_5hmC <- lapply(ls_5mC_5hmC , function(x) x[!grep("_",seqnames),]) #remove unnecasrry chroms
  
  # ls_5mC_5hmC[[grep("5mC.bed",names(ls_5mC_5hmC))]]
  
  ################
  #rename methylation data
  setnames(ls_5mC_5hmC[[grep("5mC.bed",names(ls_5mC_5hmC))]],
           "methylation_percentage", "methylation_percentage_5mC"
  )
  
  setnames(ls_5mC_5hmC[[grep("5hmC.bed",names(ls_5mC_5hmC))]],
           "methylation_percentage", "methylation_percentage_5hmC"
  )
  
  #merge 5hmc and 5mc
  Dt <- ls_5mC_5hmC[[1]][ls_5mC_5hmC[[2]], on=c("seqnames", "start", "end" ,"number_reads" ,"strand")]
  #deduce counts of methylated and unmethylated
  Dt[,`:=`(Unmethylated_percentage = (100 - (methylation_percentage_5hmC+methylation_percentage_5mC)),
           methylation_percentage_5mC_5hmC = (methylation_percentage_5hmC+methylation_percentage_5mC))][,
                                                                                                        Meth_or_UnMeth := fcase(methylation_percentage_5mC_5hmC >= methyl_percentage_threshold & number_reads >= methyl_numb_reads_threshold,"M" , default = "U")]
  #get real identity of cpg site
  Dt[, methylation_state := fcase(
    Meth_or_UnMeth == "M", 
    fifelse(methylation_percentage_5hmC > methylation_percentage_5mC, "5hmC", "5mC"),
    default = "U"
  )]
  
message("\nsanity checks\n")
print(Dt[,.N,by=methylation_state])
message("\n")
################# ##################
#make ready to export as coverage file for use in bsmooth:
# <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
  Dt[,`:=`(seqnames=seqnames, 
           start=start,end=end, 
           methylation_percentage_5mC_5hmC=methylation_percentage_5mC_5hmC, 
           number_reads=number_reads,
           count_methylated=round((number_reads*(methylation_percentage_5mC_5hmC/100))))][,count_unmethylated:=(number_reads-count_methylated)]
  return(Dt)
}
utils::globalVariables(c("seqnames", "start","end", "methylation_percentage_5mC_5hmC", 
"methylation_percentage_5hmC","methylation_percentage_5mC","lsfiles","methylation_state",
"count_methylated", "strand", "count_unmethylated", "Meth_or_UnMeth"))