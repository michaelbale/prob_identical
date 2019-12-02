# prob_identical
R Scripts for SOR index described in Patro et al., PNAS 2019

RShiny app located at https://michaelbale.shinyapps.io/prob_identical/


Pull command line version using

install.packages('devtools')
devtools::source_url("https://github.com/michaelbale/prob_identical/blob/master/prob_identical.r?raw=TRUE")


use:
SOR.obj <- getSOR(inPath, format)

Value:
S3 object containing 3 attributes:

rake_list is a list of rake IDs that show what sequences are identical to each other
p.values is a named vector in the same order as rake_list and holds the results of the test
dist.plot is a ggplot2 plot object that shows a histogram of the distance data
