#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

# load data
mRNA <- mm_mRNA

mRNA <- mRNA[1:20,]

MAE <- startObject(miR = NULL, mRNA = mRNA)

MAE <- getIdsMrnaMouse(MAE = MAE, mRNA =  assay(MAE, 2))

#test outputs have expected number of columns
test_that("mRNA_Id data have two columns", {
    expect_that(as.numeric(ncol(assay(MAE, 3))), equals(2))
    expect_that(as.numeric(ncol(assay(MAE, 4))), equals(2))
})

#internal checks
mouse <- biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl",
                             host = "useast.ensembl.org")

glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                              "entrezgene_id"),
               filters = "external_gene_name",
               values = rownames(mRNA), mart = mouse, uniqueRows = TRUE)

#check 1
#glist should be three columns
test_that("glist from biomart should be 3 columns", {
    expect_that(as.numeric(ncol(glist)), equals(3))
})

#continue
glist <- glist[! duplicated(glist$external_gene_name),]

gene_data <- cbind(mRNA, rownames(mRNA))

#check 2
#Last row should be named `rownames(mRNA)` and == to rownames(gene_data)
test_that("Final row of gene_data should be == to rownames of gene_data", {
    expect_that(rev(names(gene_data))[1], equals("rownames(mRNA)"))
    expect_equal(as.character(gene_data$`rownames(mRNA)`),
    as.character(rownames(gene_data)))
})

#contunue
#final two columns should be named 'ensembl_gene_id' and 'entrezgene_id'
m_dat <- merge(x = gene_data, y = glist, by.x = 'rownames(mRNA)',
               by.y = 'external_gene_name', all = TRUE)
test_that("last two columns in m_dat are what they are expected to be", {
    expect_that(rev(names(m_dat))[1], equals("entrezgene_id"))
    expect_that(rev(names(m_dat))[2], equals("ensembl_gene_id"))
})

#continue
rownames(m_dat) <- m_dat$`rownames(mRNA)`

MAE2 <- MultiAssayExperiment(list(
                                    mRNA_entrez = data.frame(cbind(
                                        GENENAME = rownames(m_dat),
                                        ID = m_dat$entrezgene_id)),
                                    mRNA_enembl = data.frame(cbind(
                                        GENENAME = rownames(m_dat),
                                        ID = m_dat$ensembl_gene_id))))


# MEA <- c(MAE, MAE2)
#check 3
#function output is the same as manual output
test_that("manual output == funciton output", {
    expect_identical(assay(MAE, 3), assay(MAE2, 1))
    expect_identical(assay(MAE, 4), assay(MAE2, 2))
})

#save output
saveRDS(MAE2, "IDs_mouse_mRNA.rds", compress = "xz")
