#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(clusterProfiler)
library(MultiAssayExperiment)
#download GMT mouse
#test function
MAE <- MultiAssayExperiment()
dloadGMT(MAE, speciesInitials = 'Mm') -> MAE
#internal tests
download.file("http://data.wikipathways.org/20200110/gmt/wikipathways-20200110-gmt-Mus_musculus.gmt",
"mus.gmt")
read.gmt("mus.gmt") -> musgmt
as.data.frame(musgmt) -> musgmt
colnames(musgmt) <- c("ont", "gene")
musgmt %>% tidyr::separate(ont, c("name","version","wpid",
"org"), "%") -> full_list
#check 1
#aspects of full_list
test_that("aspects of the mouse wikipathway data", {
expect_equal(length(names(full_list)), 5)
expect_equal(full_list[1,4], "Mus musculus")
})
#continue
path_gene_manual <- full_list %>% dplyr::select(wpid, gene)
path_name_manual <- full_list %>% dplyr::select(wpid, name)
#check 2
#wp entries should be the same in wpid2gene and wpid2name
test_that("wikipathways ID's are the same", {
expect_equal(path_gene_manual[,1], path_name_manual[,1])
})
#continue
path_data_manual <- full_list %>% dplyr::select(wpid, gene, name)
path_list_manual <- list(path_gene = MAE@ExperimentList$path_gene,
path_name = MAE@ExperimentList$path_name,
path_data = MAE@ExperimentList$path_data)
#check 3
#list should be 3 long
test_that("list should be 3 dataframes long", {
expect_equal(length(path_list_manual), 3)
})
#check 4
#check manual and functional outputs are the same
test_that("manual and functional output is the same", {
expect_equal(MAE@ExperimentList$path_data, path_data_manual)
expect_equal(MAE@ExperimentList$path_name, path_name_manual)
expect_equal(MAE@ExperimentList$path_gene, path_gene_manual)
})
#save data
#reduce wplists
MAE@ExperimentList$path_gene[1:2000,] -> MAE@ExperimentList$path_gene
MAE@ExperimentList$path_name[1:2000,] -> MAE@ExperimentList$path_name
MAE@ExperimentList$path_data[1:2000,] -> MAE@ExperimentList$path_data
saveRDS(MAE, "wpdata.rds", compress = "xz")
#same for human data
file.remove("mus.gmt")

