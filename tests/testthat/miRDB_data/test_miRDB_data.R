setwd("~/Documents/Package/smiRk/tests/testthat/miRDB_data/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
#load data
dloadmiRDB() -> miRDB
#check aspects of miRDB
test_that("miRDB contains miR-mRNA interaction data", {
  expect_that(length(colnames(miRDB)), equals(3))
  expect_true(is.numeric(miRDB$V3))
  expect_equal(length(list.files()), 2)
})
###############################################################################
#run function
miRDB_data(miRDB = miRDB, species = 'mmu', orgDB = org.Mm.eg.db) -> miRDB_mmu
miRDB_mmu[1:100,] -> miRDB_mmu
#check 1
#miRDB_mmu should be smaller than input and 3 columns long
test_that("miRDB_mmu is expected to be 3 columns long", {
expect_equal(length(names(miRDB_mmu)), 3)
expect_gt(length(rownames(miRDB)),
length(rownames(miRDB_mmu)))
})
#continue
#internal checks
names(miRDB) <- c('miR', 'Target', 'Score')
miRDB %>%
filter(str_detect(miR, 'mmu')) -> miRDB_s
#check 2
#miRDB_s is 3 columns long and smaller than input
test_that("miRDB_s is expected to be 3 columns long", {
expect_equal(length(names(miRDB_s)), 3)
expect_gt(length(rownames(miRDB)),
length(rownames(miRDB_s)))
})
#continue
bitr(miRDB_s$Target, fromType = 'REFSEQ', toType = 'SYMBOL',
OrgDb = org.Mm.eg.db) -> miRDB_mRNA
#check 2
#miRDB_mRNA should be 2 columns long
test_that("miRDB_mRNA is expected to be 2 columns long", {
expect_equal(length(names(miRDB_mRNA)), 2)
})
#continue
merge(miRDB_s, miRDB_mRNA, by.x = 'Target', by.y = 'REFSEQ',
all = TRUE) -> mirDB_merged
miRDB_df <- data.frame(miRDB_Interactions = paste(mirDB_merged$miR,
':', mirDB_merged$SYMBOL,sep = ''),
 miRDB_miR = mirDB_merged$miR,
 miRDB_mRNA = mirDB_merged$SYMBOL)
miRDB_df[1:100,] -> miRDB_df
#check 3
#miRDB_df should be the same as miRDB_mmu
test_that("output of function and manual code should be the same", {
expect_equal(miRDB_df, miRDB_mmu)
})
#save data
as.matrix(miRDB_df) -> miRDB_matrix
saveRDS(miRDB_matrix, "miRDB_resuts.rds", compress = "xz")

