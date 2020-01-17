#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(biomaRt)
#load data
hs_mRNA -> mRNA
mRNA[1:50,] -> mRNA
#test function
getIDs_mRNA_human(mRNA = mRNA, mirror = 'useast')
#check 1
test_that("output have two columns", {
expect_equal(length(colnames(mRNA_ensembl)), 2)
expect_equal(length(colnames(mRNA_entrez)), 2)
expect_equal(length(rownames(mRNA_entrez)),
length(rownames(mRNA_ensembl)))
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
expect_equal(length(rownames(m_dat)), 50)
expect_equal(rownames(m_dat),
m_dat$Genes)
})
#continue
as.data.frame(cbind(GENENAME = rownames(m_dat),
ID = m_dat$entrezgene_id)) -> mRNA_entrez_manual
as.data.frame(cbind(GENENAME = rownames(m_dat),
ID = m_dat$ensembl_gene_id)) -> mRNA_ensembl_manual
#check 4
test_that("functional and manual output is the same",{
expect_equal(mRNA_entrez_manual, mRNA_entrez)
expect_equal(mRNA_ensembl_manual, mRNA_ensembl)
})
#save data
saveRDS(mRNA_ensembl, "hs_mRNA_ensembl.rds", compress = "xz")
saveRDS(mRNA_entrez, "hs_mRNA_entrez.rds", compress = "xz")
