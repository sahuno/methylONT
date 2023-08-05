#' function to make aggregate statistics from per read modified probabilities from megalodon
#'
#' @param dt data.table object
#' @return aggregate statistics of per read modified bases
#' 
#' @export
#' @import data.table
#' @import utils
agg_per_read <- function(dt){
    stopifnot(nrow(dt) >= 1)
dt_out <- dt[!grep("_",chrm),][,exp_mod_log_prob:=exp(mod_log_prob)][,
                        list(PrM = sum(exp_mod_log_prob), PrUnM = 1 - sum(exp_mod_log_prob)),
                        by=list(chrm, pos, strand, read_id)][,
                        list(number_reads=.N, mean_PrUnM = mean(PrUnM), median_PrUnM = stats::median(PrUnM)),
                        by=list(chrm, pos, strand)]
return(dt_out)
}
utils::globalVariables(c("chrm", "mod_log_prob", "exp_mod_log_prob", "pos", "strand", "read_id", "PrM", "PrUnM", "number_reads", "mean_PrUnM", "median_PrUnM"))

#.datatable.aware <- TRUE

#use_package("utils")
#usethis::use_import_from("data.table", ":=")

#good resource for data.table in R packages
#https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html
# gen = function (n = 100L) {
#   dt = as.data.table(list(id = seq_len(n)))
#   dt[, grp := ((id - 1) %% 26) + 1
#      ][, grp := letters[grp]
#        ][]
# }
# aggr = function (x) {
#   stopifnot(
#     is.data.table(x),
#     "grp" %in% names(x)
#   )
#   x[, .N, by = grp]
# }
##' @examples
##' ##dt <- read.table("test.txt", header=TRUE) #read in test data 
##' agg_per_read(dt)