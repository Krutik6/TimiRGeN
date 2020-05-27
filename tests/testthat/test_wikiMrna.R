library(TimiRGeN)
library(testthat)

#load data
#negative test
log2fc <- readRDS("log2fc.rds")

path_data <- data.frame("wpid" = c(rep("WP571", 6)),
                        "gene" = c(16175, 12370,26419, 19249, 19645, 18479),
                        "name" = c(rep("Fas pathway and Stress induction of HSP regulation", 6)))

MAE <- MultiAssayExperiment()

MAE <- reduceWiki(MAE, path_data = path_data,
                  stringWiki = "Fas pathway and Stress induction of HSP regulation")

MAE <- wikiMrna(MAE,
                mRNA_express = assay(log2fc, 3),
                singleWiki = assay(MAE, 1),
                stringWiki = 'Fas pathway and Stress induction of HSP regulation')

deres <- assay(log2fc, 3)

singlewiki <- assay(MAE, 1)

GenesofInterest <- deres[which(deres$ID %in% singlewiki$gene),]

test_that("There will be no genes of interest",{
    expect_gt(1, length(rownames(GenesofInterest)))
    expect_equal(length(names(GenesofInterest)),
    length(names(deres)))
})

#positive test
test <- readRDS("interactions.rds")

GenesofInterest2 <- deres[which(deres$ID %in% test$gene),]
test_that("There will be some genes of interest",{
    expect_gt(length(rownames(GenesofInterest2)),
    length(rownames(GenesofInterest)))
    expect_equal(length(names(GenesofInterest2)),
    length(names(deres)))
})

saveRDS(GenesofInterest2, "GenesofInterest.rds", compress = "xz")
