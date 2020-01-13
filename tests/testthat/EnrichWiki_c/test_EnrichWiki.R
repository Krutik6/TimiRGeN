setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(org.Mm.eg.db)
library(clusterProfiler)
#load data
readRDS("eNames_c/elist.rds") -> elist
readRDS("downloadGMT/mouse_wplist.rds") -> mm_wp_list
mm_wp_list[[1]] -> path_gene
mm_wp_list[[2]] -> path_name
#test EnrichWiki
EnrichWiki(method = 'c', e_list = elist, orgDB = org.Mm.eg.db,
path_gene = path_gene, path_name = path_name,
ID = 'ENTREZID', universe = path_gene$gene) -> sigwiki
#internal checks
lapply(elist, function(x){
enricher(x, TERM2GENE = path_gene, TERM2NAME = path_name,
universe = path_gene$gene,
pvalueCutoff = 0.05, qvalueCutoff = 0.2,
pAdjustMethod = 'BH')}) -> lst2
#check 1
#should be a list of 5
test_that("List is of 5", {
expect_equal(length(lst2),5)
})
#continue
lapply(lst2, function(x){
setReadable(x, org.Mm.eg.db, keyType = 'ENTREZID')}) -> W_list
paste0(names(W_list), '_wikipathways', sep='') -> names(W_list)
#check 2
#check if W_list and sigwiki are the same
test_that("manual and functional output test", {
expect_equal(W_list, sigwiki)
})
#save data
saveRDS(sigwiki, "EnrichWiki_c/EnrichWiki.rds", compress = "xz")
