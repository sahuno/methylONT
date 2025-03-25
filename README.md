
<!-- README.md is generated from README.Rmd. Please edit that file -->

# methylONT

<!-- badges: start -->
<!-- badges: end -->

The goal of methylONT is to porvide utility functions to differential
methylation analysis for lonfg read sequencing

## Installation

You can install the development version of methylONT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sahuno/methylONT")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(methylONT)
#> Warning: replacing previous import 'IRanges::shift' by 'data.table::shift' when
#> loading 'methylONT'
#> Warning: replacing previous import 'Biostrings::tail' by 'utils::tail' when
#> loading 'methylONT'
#> Warning: replacing previous import 'IRanges::stack' by 'utils::stack' when
#> loading 'methylONT'
#> Warning: replacing previous import 'IRanges::relist' by 'utils::relist' when
#> loading 'methylONT'
#> Warning: replacing previous import 'Biostrings::head' by 'utils::head' when
#> loading 'methylONT'
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
