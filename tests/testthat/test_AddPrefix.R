#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load filtered_genelist
readRDS("MAE_mm.rds") -> MAE
#check 1
MAE@ExperimentList$miR -> MAE@ExperimentList$miR_p
MAE@ExperimentList$mRNA -> MAE@ExperimentList$mRNA_p

ifelse(test = grepl("miR",
names(MAE@ExperimentList$miR)) == FALSE,
yes = colnames(MAE@ExperimentList$miR) <- paste("miR",
colnames(MAE@ExperimentList$miR), sep = '_'),
no = print('miR/mRNA info is fine'))

ifelse(test = grepl("mRNA", names(
MAE@ExperimentList$mRNA)) == FALSE,
yes = colnames(MAE@ExperimentList$mRNA) <- paste("mRNA", 
colnames(MAE@ExperimentList$mRNA), sep = '_'),
no = print('miR/mRNA info is fine'))

test_that("now colnames are not equal", {
expect_false(isTRUE(all.equal(names(MAE@ExperimentList$miR),
names(MAE@ExperimentList$miR_p))))
expect_false(isTRUE(all.equal(names(MAE@ExperimentList$mRNA),
names(MAE@ExperimentList$mRNA_p))))
})

saveRDS(MAE, file = "MAE_Prefix.rds", compress = "xz")

