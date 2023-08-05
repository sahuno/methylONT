#' function to read megalodon per read modified bases and compute stats
#'
#' @param dt data.table object
#' @return aggregate statistics of per read modified bases
#' 
#' @export
#' @import data.table
#' @import utils

#to do: add script from `/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/scripts/R/aggregate_per_read_5mC_5hmc_main.r`
# function(dt){
#     message("aggregating reads")
# library(data.table)
# dt <- copy(per_mod_base_DT)
# dt2 <- dt[!grep("_",chrm),][, `:=`(exp_mod_log_prob = exp(mod_log_prob))]
# [, list(number_reads=.N,
#         mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
#         median_exp_mod_log_prob=median(exp_mod_log_prob)),by=list(chrm, pos, strand, mod_base)]

# [,`:=`(sum_Exp_mod_can_prob=exp_mod_log_prob+exp_can_log_prob, sum_mod_can_prob=mod_log_prob+can_log_prob)]
# summary(dt2$sum_mod_can_prob)
# summary(dt2$sum_Exp_mod_can_prob)                                    
# split(dt2, dt2$mod_base)

# #canonical prob + modified prob = 1 for each read
# #conclusion use the modified probs per file


#                                     # [, list(number_reads=.N,
#                                     #     mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
#                                     #     median_exp_mod_log_prob=median(exp_mod_log_prob) 
#                                     #     # min_exp_mod_log_prob=min(exp_mod_log_prob), 
#                                     #     # max_exp_mod_log_prob=max(exp_mod_log_prob),
#                                     #     # mean_exp_can_log_prob=mean(exp_can_log_prob), 
#                                     #     # median_exp_can_log_prob=median(exp_can_log_prob) 
#                                     #     # min_exp_can_log_prob=min(exp_can_log_prob), 
#                                     #     # max_exp_can_log_prob=max(exp_can_log_prob)
#                                     #     ), 
#                                     #     by=list(chrm, pos, strand, mod_base)]

# }

