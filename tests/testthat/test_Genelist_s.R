#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(stringr)
#load filtered_genelist
MAE <- readRDS(file = "MAE_Prefix.rds")
#test function
GenesList(method = "s", miR_data = MAE@ExperimentList$miR_p,
mRNA_data = MAE@ExperimentList$mRNA_p) -> MAE@metadata$genelist
#check 1
test_that("genelist should be a list of lists", {
expect_equal(class(MAE@metadata$genelist), "list")
expect_equal(class(MAE@metadata$genelist[1]), "list")
expect_equal(class(MAE@metadata$genelist[2]), "list")
expect_equal(class(MAE@metadata$genelist[1][1]), "list")
})
#internal checks
genedata <- list(miR_data = MAE@ExperimentList$miR_p,
mRNA_data = MAE@ExperimentList$mRNA_p)
#check 2
test_that("genedata is a list", {
expect_equal(class(genedata), "list")
})
#contineu
genes.split <- lapply(genedata, function(df)
{
Unigenes <- unique(str_extract(names(df),"\\S+\\."))
List <- lapply(Unigenes,function(name){return(df[,grep(name,names(df),
fixed=TRUE)])})
names(List) <- Unigenes
return(List)
})
#check 3
test_that("internal and external output is the same", {
expect_equal(genes.split, MAE@metadata$genelist)
})
#save data
saveRDS(MAE@metadata$genelist, "genelist_s.rds")

