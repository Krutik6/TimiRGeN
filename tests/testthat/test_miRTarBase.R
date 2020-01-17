#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(tidyverse)
#load data
dloadmiRTarBase() -> miRTarBase
#checks
#check aspects of miRTarBase
test_that("miRTarBase contains miR-mRNA interaction data,
would expect fewer strongly assesed interactions and more
weakly assessed interactions", {
expect_true(is.data.frame(miRTarBase))
expect_equal(names(miRTarBase)[2], "miRNA")
expect_equal(names(miRTarBase)[4], "Target.Gene")
expect_lt(length(rownames(miRTarBase[which(miRTarBase$Support.Type
== 'Functional MTI'),])),
length(rownames(miRTarBase[which(miRTarBase$Support.Type
== 'Functional MTI (Weak)')
,])))
})
#remove miRTarBase.csv
file.remove("miRTarBase.csv")
#run function
miRTarBase_data(mirtarbase = miRTarBase, species = 'mmu') -> miRTarBase_mmu
miRTarBase_mmu[1:100,] -> miRTarBase_mmu
#check 1
#miRTarBase_mmu should be smaller than input and 3 columns long
test_that("miRTarBase_mmu is expected to be 3 columns long", {
expect_equal(length(names(miRTarBase_mmu)), 3)
expect_gt(length(rownames(miRTarBase)),
length(rownames(miRTarBase_mmu)))
})
#internal checks
miRTarBase %>%
filter(str_detect(miRNA, 'mmu')) -> miRTarBase_s
miRTarBase_s[which(miRTarBase_s$Support.Type == 'Functional MTI')
,] -> miRTarBase_Fun
#check 2
#should be fewer in mIRTarBase_s than in input and even less in miRTarBase_Fun
test_that("Progresively ferwer rows", {
expect_lt(length(rownames(miRTarBase_s)),
length(rownames(miRTarBase)))
expect_lt(length(rownames(miRTarBase_Fun)),
length(rownames(miRTarBase_s)))
})
#continue
miRTarBase_df <- data.frame(miRTarBase_Interactions = paste(
miRTarBase_Fun$miRNA,':', miRTarBase_Fun$Target.Gene, sep = ''),
miRTarBase_microRNA = miRTarBase_Fun$miRNA,
miRTarBase_mRNA = miRTarBase_Fun$Target.Gene)
miRTarBase_df[1:100,] -> miRTarBase_df
#check 3
#output should be the same as function output
test_that("output should be the same", {
expect_equal(miRTarBase_df, miRTarBase_mmu)
})
#save data
as.matrix(miRTarBase_df) -> miRTarBase_matrix
saveRDS(miRTarBase_matrix, "miRTarBase_results.rds", compress = "xz")
