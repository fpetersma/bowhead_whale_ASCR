## =============================================================================
## SCRIPT TO UPDATE PACKAGE FOR RCPP FUNCTIONALITY
## =============================================================================

library(ascrRcpp)
detach("package:ascrRcpp", unload=TRUE)

library(Rcpp)

compileAttributes("Scripts/bowhead whales/ascrRcpp")

remove.packages("ascrRcpp")

# devtools::build("../ascrRcpp")
devtools::install("Scripts/bowhead whales/ascrRcpp")
