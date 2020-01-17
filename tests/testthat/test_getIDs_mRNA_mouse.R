#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(biomaRt)
# load data
mm_mRNA -> mRNA
mRNA[1:50,] -> mRNA
getIDs_mRNA_mouse(mRNA)
#test outputs have expected number of columns
test_that("mRNA_Id data have two columns", {
expect_that(as.numeric(ncol(mRNA_entrez)), equals(2))
expect_that(as.numeric(ncol(mRNA_ensembl)), equals(2))
})
#internal checks
mouse <- biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl",
GRCh=37, host = "useast.ensembl.org")
glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
"entrezgene_id"), filters = "external_gene_name",
values = rownames(mRNA), mart = mouse, uniqueRows = TRUE)
#check 1
#glist should be three columns
test_that("glist from biomart should be 3 columns", {
expect_that(as.numeric(ncol(glist)), equals(3))
})
#continue
glist[! duplicated(glist$external_gene_name),] -> glist
cbind(mRNA, rownames(mRNA)) -> gene_data
#check 2
#Last row should be named `rownames(mRNA)` and == to rownames(gene_data)
test_that("Final row of gene_data should be == to rownames of gene_data", {
expect_that(rev(names(gene_data))[1], equals("rownames(mRNA)"))
expect_equal(as.character(gene_data$`rownames(mRNA)`),
as.character(rownames(gene_data)))
})
#contunue
#final two columns should be named 'ensembl_gene_id' and 'entrezgene_id'
merge(x = gene_data, y = glist, by.x = 'rownames(mRNA)',
by.y = 'external_gene_name', all = TRUE) -> m_dat
test_that("last two columns in m_dat are what they are expected to be", {
expect_that(rev(names(m_dat))[1], equals("entrezgene_id"))
expect_that(rev(names(m_dat))[2], equals("ensembl_gene_id"))
})
#continue
rownames(m_dat) <- m_dat$`rownames(mRNA)`
as.data.frame(cbind(GENENAME = rownames(m_dat),
ID = m_dat$entrezgene_id))-> mRNA_entrez_man
as.data.frame(cbind(GENENAME = rownames(m_dat),
ID = m_dat$ensembl_gene_id)) -> mRNA_ensembl_man
#check 3
#function output is the same as manual output
test_that("manual output == funciton output", {
expect_identical(mRNA_ensembl, mRNA_ensembl_man)
expect_identical(mRNA_entrez, mRNA_entrez_man)
})

#save output
saveRDS(mRNA_ensembl, "mm_mRNA_ensembl.rds")
saveRDS(mRNA_entrez, "mm_mRNA_entrez.rds")



