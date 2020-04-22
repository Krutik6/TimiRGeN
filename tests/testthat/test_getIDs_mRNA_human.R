#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
library(biomaRt)
#load data
mRNA <- hs_mRNA
mRNA <- mRNA[1:20,]
MAE <- startObject(miR = NULL, mRNA = mRNA)
#test function
MAE <- getIdsMrnaHuman(MAE, mRNA = assay(MAE, 2),
                       mirror = 'useast')
#check 1
test_that("output have two columns", {
    expect_equal(length(colnames(assay(MAE, 3))), 2)
    expect_equal(length(colnames(assay(MAE, 4))), 2)
    expect_equal(length(rownames(assay(MAE, 3))),
    length(rownames(assay(MAE, 4))))
})
#internal checks
mRNA$Genes <- rownames(mRNA)
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
m_dat <- merge(x = mRNA, y = glist, by.x = 'Genes',
               by.y = 'external_gene_name', all = TRUE)

m_dat <- m_dat[! duplicated(m_dat$Genes),]
m_dat <- m_dat[order(m_dat$Genes),]
rownames(m_dat) <- m_dat$Genes
#check 3
test_that("Genes and rownames are the same", {
    expect_equal(length(colnames(m_dat)), 9)
    expect_equal(length(rownames(m_dat)), 20)
    expect_equal(rownames(m_dat),m_dat$Genes)
})
#continue
MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                                mRNA_entrez = data.frame(cbind(
                                                    GENENAME = rownames(m_dat),
                                                    ID = m_dat$entrezgene_id)),
                                                mRNA_enembl = data.frame(cbind(
                                                    GENENAME = rownames(m_dat),
                                                    ID = m_dat$ensembl_gene_id)
    ))))
#check 4
test_that("functional and manual output is the same",{
    expect_equal(assay(MAE, 3), assay(MAE2, 1))
    expect_equal(assay(MAE, 4), assay(MAE2, 2))
})
