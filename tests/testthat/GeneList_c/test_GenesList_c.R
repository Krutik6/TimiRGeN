setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load genetic_data
readRDS("CombineGenes/genetic_data.rds") -> genetic_data
#perform GenesList combined function
GenesList(method = 'c', genetic_data = genetic_data,
timeString = 'D') -> geneslist
#internal checks
gsub(x = colnames(genetic_data), pattern = 'D',
replacement = 'TP') -> colnames(genetic_data)
#check1
#genetic_data should not have 'D' in the colnames
test_that("genetic_data colnames has not 'D' in it", {
for (i in 1:10) {
expect_false(grepl(colnames(genetic_data), pattern = 'D')[i])
}
})
#continue
lapply(split.default(genetic_data, sub("(TP\\d+).*", "\\1",
names(genetic_data))), as.list) -> X
#check 2
#should be a list of 5 lists
test_that("X is a list of 5", {
expect_equal(length(X), 5)
})
#continue
gsub(names(X), pattern = 'TP', replacement = 'D') -> names(X)
#check 3
#there should be 'D' in list names
test_that("genetic_data colnames has not 'D' in it", {
for (i in 1:5) {
expect_true(grepl(names(X), pattern = 'D')[i])
}
})
#continue
L1 <- lapply(X, data.frame, stringsAsFactors = FALSE)
L2 <- lapply(L1, function(DF) {rownames(DF) <- rownames(genetic_data); DF})
#check 4
#nested DFs in L2 have the same rownames as genetic_data
test_that("Dfs in L2 have the same rownames as in genetic_data", {
lapply(L2, function(x) { rownames(x)}) -> a
for (i in 1:5) {
expect_equal(as.character(unlist(a[i])),
as.character(rownames(genetic_data)))
}
})
#continue
gtools::mixedsort(names(L2)) -> C1
L2 <- L2[C1]
L3 <- L2[gtools::mixedsort(names(L2))]
genelist_man <- lapply(L3, function(x) {
colnames(x) <- sub(x = colnames(x), pattern = 'TP', replacement = 'D')
x}
)
#check 5
#are genelist_man and genelist the same?
test_that("manual and function output is the same", {
expect_equal(genelist_man, geneslist)
})
#savedata
saveRDS(geneslist,"GeneList_c/geneslist.rds", compress = "xz")
