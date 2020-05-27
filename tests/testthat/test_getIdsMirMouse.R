#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(org.Mm.eg.db)

# load data
miR <- mm_miR

miR <- miR[1:20,]

MAE <- startObject(miR = miR, mRNA = NULL)

#test function
MAE <- getIdsMirMouse(MAE, assay(MAE, 1))

#test outputs have expected number of columns
test_that("miR_Id data have two columns", {
    expect_that(as.numeric(ncol(assay(MAE, 3))), equals(2))
    expect_that(as.numeric(ncol(assay(MAE, 4))), equals(2))
    expect_that(as.numeric(ncol(assay(MAE, 5))), equals(2))
    expect_that(as.numeric(ncol(assay(MAE, 6))), equals(2))
})

#internal checks
miR$Genes <- miR$names <- rownames(miR)

miR$Genes <- sub(x = miR$Genes, pattern = "-3p", replacement = "")

miR$Genes <- sub(x = miR$Genes, pattern = "-5p", replacement = "")

#check 1
#test that -3p and -5p have been removed
test_that("-3p and -5p should be absent", {
    expect_that(length(grep(miR$Genes, pattern = "-5p")), equals(0))
    expect_that(length(grep(miR$Genes, pattern = "-3p")), equals(0))
})

#continue
miR$Genes <- micrornaFull(miRdf = miR$Genes, species = 'mmu')

#check 2
test_that("Genes and names are different", {
    expect_false(isTRUE(all.equal(miR$Genes, miR$names)))
})

#continue
miR_entrez2 <- bitr(geneID = miR$Genes, fromType = 'GENENAME',
                    toType = 'ENTREZID',OrgDb = org.Mm.eg.db)

miR_ensembl2 <- bitr(geneID = miR$Genes, fromType = 'GENENAME',
                     toType = 'ENSEMBL', OrgDb = org.Mm.eg.db)

#check 3
test_that("miR_entrez and miR_ensembl are different but similiar", {
    expect_equal(length(names(miR_ensembl2)), 2)
    expect_equal(length(names(miR_entrez2)), 2)
})

#conitnue
miR_merged <- merge(x = miR, y = miR_ensembl2, by.x = 'Genes',
                    by.y = 'GENENAME', all = TRUE)

miR_merged <- merge(x = miR_merged, y = miR_entrez2, by.x = 'Genes',
                    by.y = 'GENENAME', all = TRUE)

miR_merged <- miR_merged[!duplicated(miR_merged$names),]

miR_merged <- miR_merged[order(miR_merged$names),]

#check 4
test_that("miR_merged has 14 columns", {
    expect_equal(length(names(miR_merged)), 14)
})

#continue
miR_merged$ENTREZID_adjusted <- nonUnique(Col = miR_merged$ENTREZID,
                                           sep = ".",suffix = "")

miR_merged$ENSEMBL_adjusted <- nonUnique(Col = miR_merged$ENSEMBL,
                                          sep = ".", suffix = "")

#check 4
test_that("miR_merged has 16 columns", {
    expect_equal(length(names(miR_merged)), 16)
})

rownames(miR_merged) <- miR_merged$names

#continue
MAE2 <- MultiAssayExperiment(list(
    miR_entrez = data.frame(cbind(
                            GENENAME = rownames(miR_merged),
                            ID = miR_merged$ENTREZID)),
    miR_ensembl = data.frame(cbind(
                            GENENAME = rownames(miR_merged),
                            ID = miR_merged$ENSEMBL)),
    miR_adjusted_entrez = data.frame(cbind(
                            GENENAME = rownames(miR_merged),
                            ID = miR_merged$ENTREZID_adjusted)),
    miR_adjusted_ensembl = data.frame(cbind(
                            GENENAME = rownames(miR_merged),
                            ID = miR_merged$ENSEMBL_adjusted))))


# MAE <- c(MAE, MAE2)
test_that("manual and functional output are the same", {
    expect_equal(assay(MAE, 3), assay(MAE2, 1))
    expect_equal(assay(MAE, 4), assay(MAE2, 2))
    expect_equal(assay(MAE, 5), assay(MAE2, 3))
    expect_equal(assay(MAE, 6), assay(MAE2, 4))

})

#save output
saveRDS(MAE2, "IDs_mouse_miR.rds", compress = "xz")
