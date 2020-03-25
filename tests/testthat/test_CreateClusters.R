#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(Mfuzz)
library(MultiAssayExperiment)
library(rWikiPathways)
#load data
Pmat <- as.data.frame(readRDS("Pmat.rds"))
MAE <- MultiAssayExperiment(list("Pmat" = Pmat))

#test function
MAE <- CreateClusters(MAE, method = 'c', 
                      percentMatrix = as.matrix(assay(MAE, 1)),
                      noClusters = 5, variance = 0)


dev.off()
#test output
#check 1
test_that("the output is in the expected form", {
    expect_equal(length(colnames(assay(MAE, 1))), 5)
    expect_true(is.list(metadata(MAE)))
})
#internal checks
df <- as.data.frame(t(assay(MAE, 1)))
#check 2
test_that("df has the correct number of columns", {
    expect_equal(length(rownames(assay(MAE, 1))), length(colnames(df)))
})
#continue
df$Total <- NULL
df[vapply(df, is.factor, logical(1))] <- lapply(df[vapply(df, is.factor,
               logical(1))], function(x) as.numeric(as.character(x)))

df <- round(df, 0)
df <- na.omit(df)
df2 <- df
Eset <- new('ExpressionSet', exprs = as.matrix(df2))
Eset_sd <- filter.std(Eset, min.std = 0.5)
dev.off()
Eset_st <- standardise(Eset_sd)
m <- mestimate(Eset_st)
cl <- mfuzz(Eset_st, centers = 5, m=m)
X <- as.data.frame(cl$membership)
X <- head(X)
 for (i in seq_along(rownames(X))) {
    getPathwayInfo(rownames(X)[i])[[3]]
    } -> X$Description[i]
#check 3
test_that("Manual and code output is the same", {
        expect_equal(length(colnames(X)),
        length(colnames(assay(MAE, 2))))
})

#save data
saveRDS(MAE, "Clusters.rds", compress = "xz")