setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load filtered_genelist
readRDS("non_unique/miR_IDs.rds") -> miR_IDs

miR_IDs -> Y2

non_unique(Col = miR_IDs$Hs_n, sep = '-', suffix = 'p') -> miR_IDs$Hs_n

expect_false(isTRUE(all.equal(Y2$Hs_n, miR_IDs$Hs_n)))
