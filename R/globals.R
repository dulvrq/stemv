
# to prevent Notes in R CMD check
# see https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887

utils::globalVariables(c("region", "species", ".", "name", "D_lower"))
