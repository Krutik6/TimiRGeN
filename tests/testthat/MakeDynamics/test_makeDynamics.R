setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("Express/miR_log2fc.rds") -> miR_express
readRDS("Express/mRNA_log2fc.rds") -> mRNA_express
readRDS("getIDs_miR_mouse/miR_adjusted_entrez.rds") -> miR_adjusted_entrez
#test MakeDynamics
MakeDynamic(miR_expression = miR_express,
mRNA_expression = mRNA_express,
miR_IDs_adj = miR_adjusted_entrez,
Datatype = 'L') -> Dynamics
#check 1
#Dynamics should have 7 columns
test_that("Dynamics has the expected features", {
expect_equal(length(names(Dynamics)), 7)
expect_equal(length(rownames(Dynamics)),
length(rownames(miR_express)) +
length(rownames(mRNA_express)))
})
#save data
saveRDS(Dynamics, "MakeDynamics/Dynamics.rds", compress = "xz")
