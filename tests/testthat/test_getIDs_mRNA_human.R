#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
library(biomaRt)
#load data
hs_mRNA -> mRNA
mRNA[1:20,] -> mRNA
StartObject(miR = NULL, mRNA = mRNA) -> MAE
#test function
getIDs_mRNA_human(MAE, mRNA = MAE@ExperimentList$mRNA, 
mirror = 'useast') -> MAE
#check 1
test_that("output have two columns", {
expect_equal(length(colnames(MAE@ExperimentList$mRNA_ensembl)), 2)
expect_equal(length(colnames(MAE@ExperimentList$mRNA_entrez)), 2)
expect_equal(length(rownames(MAE@ExperimentList$mRNA_entrez)),
length(rownames(MAE@ExperimentList$mRNA_ensembl)))
})
#internal checks
rownames(mRNA) -> mRNA$Genes
human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
GRCh=37, host = "useast.ensembl.org")
glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
"entrezgene_id"),
filters = "external_gene_name", values = mRNA$Genes,
mart = human, uniqueRows = TRUE)
#check 2
test_that("glist has three columns", {
expect_equal(length(names(glist)), 3)
})
#continue
merge(x = mRNA, y = glist, by.x = 'Genes',
by.y = 'external_gene_name', all = TRUE) -> m_dat
m_dat[! duplicated(m_dat$Genes),] -> m_dat
m_dat[order(m_dat$Genes),] -> m_dat
rownames(m_dat) <- m_dat$Genes
#check 3
test_that("Genes and rownames are the same", {
expect_equal(length(colnames(m_dat)), 9)
expect_equal(length(rownames(m_dat)), 20)
expect_equal(rownames(m_dat),
m_dat$Genes)
})
#continue
MAE2 <- MultiAssayExperiment()
MAE2@ExperimentList$mRNA_entrez <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$entrezgene_id))
MAE2@ExperimentList$mRNA_ensembl <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$ensembl_gene_id))
#check 4
test_that("functional and manual output is the same",{
expect_equal(MAE2@ExperimentList$mRNA_entrez, MAE@ExperimentList$mRNA_entrez)
expect_equal(MAE2@ExperimentList$mRNA_ensembl, MAE@ExperimentList$mRNA_ensembl)
})
#save data
saveRDS(MAE, "IDs_human_mRNA.rds", compress = "xz")

