% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cpgs.R
\name{cpg_sites}
\alias{cpg_sites}
\title{function to make aggregate statistics from per read modified probabilities from megalodon}
\usage{
cpg_sites(species = "human", genome = "hg38", offset = 0)
}
\arguments{
\item{species}{which species to use, options "human" or "mouse"; implement for mouse}

\item{genome}{which genome to use, options "hg19" or "hg38"}

\item{offset}{0/1-based granges object}
}
\value{
granges object of cpg loci
}
\description{
function to make aggregate statistics from per read modified probabilities from megalodon
}
\examples{
cpg_sites(species="human",genome="hg38", offset=0)
}
