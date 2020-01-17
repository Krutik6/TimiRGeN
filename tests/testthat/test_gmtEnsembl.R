#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(org.Mm.eg.db)
library(clusterProfiler)
#load wpid2pathway data
readRDS("mouse_wplist.rds") -> mm_wpid2pathway
mm_wpid2pathway[[1]] -> path_gene
mm_wpid2pathway[[3]] -> path_data
#test GMT_ensembl function
GMT_ensembl(path_gene = path_gene, path_data = path_data, orgDB = org.Mm.eg.db)
#internal checks
bitr(geneID = path_gene$gene, fromType = 'ENTREZID', toType = 'ENSEMBL',
     OrgDb = org.Mm.eg.db) -> ensembl_values
#check 1
#ensembl values and entrezIDs
test_that("ensembl_values has 2 columns", {
  expect_that(length(names(ensembl_values)), equals(2))
})
#continue
merge(x = path_gene, y = ensembl_values, by.x = 'gene',
      by.y = 'ENTREZID') -> path_gene_ensembl_manual
#check 2
#wpid2gene_ensembl_manual have three columns
test_that("wpid2gene_ensembl_manual has 2 columns", {
  expect_that(length(names(path_gene_ensembl_manual)), equals(3))
})
#continue
merge(x = path_gene_ensembl_manual, y = path_data, by.x = 'gene',
      by.y = 'gene') -> path_data_ensembl_manual
#check 3
#wp_data_ensembl_manual have five columns
test_that("wp_data_ensembl_manual has 5 columns", {
  expect_that(length(names(path_data_ensembl_manual)), equals(5))
})
#continue
path_data_ensembl_manual[path_data_ensembl_manual[
  ,2] == path_data_ensembl_manual[,4],] -> path_data_ensembl_corrected
#check 4
#column 2 and 4 are the same
test_that("Columns 2 and 4 of DFs are the same but were different", {
  expect_equal(path_data_ensembl_corrected[,2], path_data_ensembl_corrected[,4])
  expect_false(isTRUE(path_data_ensembl_manual[,2]),
               equals(path_data_ensembl_manual[,4]))
})
#continue
path_gene_ensembl_manual$gene <- NULL
path_data_ensembl_corrected$gene <- path_data_ensembl_corrected$wpid.y <- NULL
names(path_gene_ensembl_manual)[2] <- 'gene'
names(path_data_ensembl_corrected)[2] <- 'gene'
#check 5
#colnames are the same between dataframes
test_that("Column features are as expected", {
  expect_equal(names(path_data_ensembl_corrected), names(path_data_ensembl))
  expect_equal(names(path_gene_ensembl_manual), names(path_gene_ensembl))
  expect_equal(length(names(path_data_ensembl_corrected)), 3)
  expect_equal(length(names(path_gene_ensembl_manual)), 2)
})
#check 6
#There are no duplicates in wp_data_ensembl_manual
path_data_ensembl_corrected[apply(path_data_ensembl_corrected, 1,
                      function(x) length(unique(x))) > 1,] -> X
test_that("X and wp_data_ensembl_manual are the same", {
  expect_equal(X, path_data_ensembl_corrected)
})
#check 7
#test manual and functional outputs
test_that("outputs are the same in function and manual", {
  expect_equal(path_gene_ensembl_manual, path_gene_ensembl)
  expect_equal(path_data_ensembl_corrected, path_data_ensembl)
})
#save file
list(path_gene_ensembl = path_gene_ensembl,
     path_data_ensembl = path_data_ensembl) -> wpath_mouse_en
saveRDS(wpath_mouse_en, "wpath_mouse_en.rds")

