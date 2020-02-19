#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(biomaRt)
library(MultiAssayExperiment)
# load data
mm_mRNA -> mRNA
mRNA[1:20,] -> mRNA
StartObject(miR = NULL, mRNA = mRNA) -> MAE
getIDs_mRNA_mouse(MAE = MAE, mRNA =  MAE@ExperimentList$mRNA) -> MAE
#test outputs have expected number of columns
test_that("mRNA_Id data have two columns", {
expect_that(as.numeric(ncol(MAE@ExperimentList$mRNA_entrez)), equals(2))
expect_that(as.numeric(ncol(MAE@ExperimentList$mRNA_ensembl)), equals(2))
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
MAE2 <- MultiAssayExperiment()
MAE2@ExperimentList$mRNA_entrez <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$entrezgene_id)) 
MAE2@ExperimentList$mRNA_ensembl <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$ensembl_gene_id))
#check 3
#function output is the same as manual output
test_that("manual output == funciton output", {
expect_identical(MAE@ExperimentList$mRNA_ensembl, MAE2@ExperimentList$mRNA_ensembl)
expect_identical(MAE@ExperimentList$mRNA_entrez, MAE2@ExperimentList$mRNA_entrez)
})

#save output
saveRDS(MAE, "IDs_mouse_mRNA.rds", compress = "xz")




