#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(org.Hs.eg.db)
#load data
miR <- hs_miR
miR <- miR[1:20,]
#keep in vignette to show that the user may have to alter their input files
rownames(miR) <- gsub(rownames(miR), pattern = "\\.", replacement =  "-")
MAE <- startObject(miR = miR, mRNA = NULL)
#test function
MAE <- getIdsMirHuman(MAE, MAE@ExperimentList$miR)
#check 1
test_that("output files should have two columns and the same number of
rownames", {
    expect_equal(length(names(assay(MAE, 3))), 2)
    expect_equal(length(names(assay(MAE, 4))), 2)
    expect_equal(length(names(assay(MAE, 5))), 2)
    expect_equal(length(names(assay(MAE, 6))), 2)
    expect_equal(names(assay(MAE, 3)), names(assay(MAE, 5)))
    expect_equal(length(rownames(assay(MAE, 4))),
    length(rownames(assay(MAE, 6))))
})
#internal checks
miR$Genes <- miR$MicroRNA <- rownames(miR)
gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "") -> miR$MicroRNA
gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "") -> miR$MicroRNA
#check 2
test_that("miR has 8 columns and the Genes and MicroRNA columns are different",
{
    expect_equal(length(names(miR)), 8)
})

#continue
miR$MicroRNA <- micrornaFull(miRdf = miR$MicroRNA,
                              species = 'hsa')

miR_entrez_manual <- bitr(geneID = miR$MicroRNA, fromType = 'GENENAME',
                          toType = 'ENTREZID',OrgDb = org.Hs.eg.db)

miR_ensembl_manual <- bitr(geneID = miR$MicroRNA, fromType = 'GENENAME',
                           toType = 'ENSEMBL',OrgDb = org.Hs.eg.db)

miR_merged <- merge(x = miR, y = miR_ensembl_manual, by.x = 'MicroRNA',
                    by.y = 'GENENAME', all = TRUE)

miR_merged <- merge(x = miR_merged, y = miR_entrez_manual,
                    by.x = 'MicroRNA', by.y = 'GENENAME', all = TRUE)

miR_merged <- miR_merged[!duplicated(miR_merged$Genes),]
miR_merged <- miR_merged[order(miR_merged$Genes),]
#check 3
test_that("miR_merged has similiar entries to miR", {
    expect_equal(length(rownames(miR_merged)),
    length(rownames(miR)))
    expect_equal( miR_merged$Genes, rownames(miR))
})
#continue
miR_merged$ENTREZID_adjusted <- nonUnique(Col = miR_merged$ENTREZID,
                                           sep = ".", suffix = "")

miR_merged$ENSEMBL_adjusted <- nonUnique(Col = miR_merged$ENSEMBL,
                                          sep = ".", suffix = "")

miR_merged <- miR_merged[! duplicated(miR_merged$Genes),]

miR_merged <- miR_merged[order(miR_merged$Genes),]
rownames(miR_merged)  <- miR_merged$Genes
#check 4
test_that("miR_merged has adjusted information", {
    expect_equal(length(names(miR_merged)), 12)
    expect_equal(names(miR_merged)[11], "ENTREZID_adjusted")
    expect_equal(names(miR_merged)[12], "ENSEMBL_adjusted")
})
#continue


MAE2 <- suppressMessages(MultiAssayExperiment(list(
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
                                        ID = miR_merged$ENSEMBL_adjusted))
)))
#check 5
test_that("manual and functional output are the same", {
    expect_equal(assay(MAE, 3), assay(MAE2, 1))
    expect_equal(assay(MAE, 4), assay(MAE2, 2))
    expect_equal(assay(MAE, 5), assay(MAE2, 3))
    expect_equal(assay(MAE, 6), assay(MAE2, 4))
})
