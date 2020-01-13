setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(biomaRt)
#load data
mm_mRNA -> mRNA
mRNA[1:50,] -> mRNA
#test function
getIDs_mRNA_mousetohuman(mRNA)
#check 1
test_that("mRNA output is as expected", {
expect_equal(length(names(mRNA_ensembl)), 2)
expect_equal(length(names(mRNA_entrez)), 2)
expect_equal(length(names(mRNA_human_renamed)), 10)
})
#internal checks
musGenes <- rownames(mRNA)
human = useEnsembl(biomart = "ensembl",
dataset = "hsapiens_gene_ensembl", mirror = 'useast')
mouse = useEnsembl(biomart = "ensembl",
dataset = "mmusculus_gene_ensembl", mirror = 'useast')
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
values = musGenes ,
mart = mouse, attributesL = c("hgnc_symbol",
"entrezgene_id",
"ensembl_gene_id"),
martL = human, uniqueRows=TRUE)
#check 2
test_that("biomart output has four columns", {
expect_that(length(names(genesV2)), equals(4))
})
#continue
genesV2[! duplicated(genesV2$MGI.symbol),] -> mh
mh[! duplicated(mh$HGNC.symbol),] -> mhu
mhu[order(mhu$MGI.symbol),] -> mouse_human
mouse_human[which(mouse_human$MGI.symbol %in% rownames(mRNA) == TRUE),] -> mh_f
#check 3
test_that("dataframes have progressively fewer rows", {
expect_gt(length(rownames(genesV2)),
length(rownames(mh)))
})
#continue
mRNA[which(rownames(mRNA) %in% mh_f$MGI.symbol == TRUE),] -> mRNA_f
mRNA_f[order(rownames(mRNA_f)),] -> mRNA_o
#check 4
test_that("mRNA_o and mh_f have the same columns",{
expect_equal(rownames(mRNA_f), mh_f$MGI.symbol)
})
#continue
cbind(mRNA_o, mh_f) -> mRNA_data
rownames(mRNA_data) <- mRNA_data$HGNC.symbol
as.data.frame(cbind(GENENAME = rownames(mRNA_data),
ID = mRNA_data$NCBI.gene.ID)) -> mRNA_entrez_manual
as.data.frame(cbind(GENENAME = rownames(mRNA_data),
ID = mRNA_data$Gene.stable.ID)) -> mRNA_ensembl_manual
mRNA_data$MGI.symbol <- mRNA_data$HGNC.symbol <- mRNA_data$NCBI.gene.ID <-
mRNA_data$Gene.stable.ID <- NULL
#check 5
test_that("manual and funcitonal outputs are the same", {
expect_equal(dim(mRNA_ensembl), dim(mRNA_ensembl_manual))
expect_equal(dim(mRNA_entrez), dim(mRNA_entrez_manual))
expect_equal(dim(mRNA_human_renamed), dim(mRNA_data))
})

