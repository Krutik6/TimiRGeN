#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
readRDS("log2fc.rds") -> Log2FC
readRDS("IDs_mouse_miR.rds") -> miR
#test MakeDynamics
MakeDynamic(miR_expression = Log2FC@ExperimentList$miR_log2fc,
mRNA_expression = Log2FC@ExperimentList$mRNA_log2fc,
miR_IDs_adj = miR@ExperimentList$miR_adjusted_entrez,
Datatype = 'L') -> Dynamics
# interal checks
Log2FC@ExperimentList$miR_log2fc -> miR_log2fc
miR@ExperimentList$miR_adjusted_entrez -> miR_adjusted_entrez
Log2FC@ExperimentList$mRNA_log2fc -> mRNA_log2fc
rownames(miR_log2fc) -> miR_log2fc$names
merge(x= miR_log2fc,
y= miR_adjusted_entrez,
by.x= 'names',
by.y= 'GENENAME', all = TRUE) -> X
rownames(X) <- X$names
X$names <- X$ID.x <- NULL
gsub(colnames(X), pattern = "ID.y", replacement = "ID") -> colnames(X)
names(Log2FC@ExperimentList$mRNA_log2fc) <- names(X)
rbind(X, Log2FC@ExperimentList$mRNA_log2fc) -> Dynamic
cbind(Dynamic, 'L') -> Dynamics2
#check 1
#Dynamics should have 7 columns
test_that("Dynamics has the expected features", {
expect_equal(length(names(Dynamics)), 7)
expect_equal(length(rownames(Dynamics)),
length(rownames(Log2FC@ExperimentList$mRNA_log2fc)) +
length(rownames(Log2FC@ExperimentList$miR_log2fc)))
})

