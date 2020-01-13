setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(org.Mm.eg.db)
library(clusterProfiler)
#load data
readRDS("eNames_s/elist_s.rds") -> elist
readRDS("downloadGMT/mouse_wplist.rds") -> mm_wp_list
mm_wp_list[[1]] -> path_gene
mm_wp_list[[2]] -> path_name
mm_wp_list[[3]] -> path_data
#GMT_ensembl
GMT_ensembl(path_gene = path_gene, path_data = path_data,
orgDB = org.Mm.eg.db)

elist[6:10] -> elist2


#test function
EnrichWiki(method = 's', e_list = elist2, orgDB = org.Mm.eg.db,
path_gene = path_gene_ensembl, path_name = path_name, ID = 'ENSEMBL',
universe = path_gene_ensembl$gene) -> sigwiki

#check 1
test_that("sigwiki has 5 lists", {
expect_equal(class(sigwiki), "list")
expect_equal(length(names(sigwiki)), 5)
})
