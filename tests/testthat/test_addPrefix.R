#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load filtered_genelist
MAE <- readRDS("MAE_mm.rds")
#check 1
miR_p <- assay(MAE, 1)

ifelse(test = grepl("miR", names(miR_p)) == FALSE,
                yes = colnames(miR_p) <- paste("miR",
                colnames(miR_p), sep = '_'),
                no = print('miR/mRNA info is fine'))

mRNA_p <- assay(MAE, 2)

ifelse(test = grepl("mRNA", names(mRNA_p)) == FALSE,
                yes = colnames(mRNA_p) <- paste("mRNA",
                colnames(mRNA_p), sep = '_'),
                no = print('miR/mRNA info is fine'))

MAE2 <- startObject(miR = miR_p, mRNA = mRNA_p)
names(MAE2) <- c("miR_p", "mRNA_p")
MAE <- c(MAE, MAE2)

test_that("now colnames are not equal", {
                expect_false(isTRUE(all.equal(names(assay(MAE, 1)),
                                names(assay(MAE, 3)))))
                expect_false(isTRUE(all.equal(names(assay(MAE, 2)),
                                names(assay(MAE, 4)))))
})

prefixString <- "mRNA"
paste(prefixString, "p", sep = "_") -> x
x[1]

MAE3 <- MultiAssayExperiment(experiments = list(x = mRNA_p))
names(MAE3) <- x
saveRDS(MAE, file = "MAE_Prefix.rds", compress = "xz")
