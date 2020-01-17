#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load filtered_genelist
mm_miR -> miR
miR[1:100,] -> miR
mm_mRNA -> mRNA
mRNA[1:200,] -> mRNA
#check 1
miR -> miR1
mRNA -> mRNA1
ifelse(test = grepl("miR", names(miR)) == FALSE,
yes = colnames(miR) <- paste("miR", colnames(miR),
sep = '_'),
no = print('miR/mRNA info is fine'))

ifelse(test = grepl("mRNA", names(mRNA)) == FALSE,
yes = colnames(mRNA) <- paste("mRNA", colnames(mRNA),
sep = '_'),
no = print('miR/mRNA info is fine'))
test_that("now colnames are not equal", {
expect_false(isTRUE(all.equal(names(miR1), names(miR))))
expect_false(isTRUE(all.equal(names(mRNA1), names(mRNA))))
})


