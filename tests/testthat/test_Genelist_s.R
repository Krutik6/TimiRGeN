#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(stringr)
#load filtered_genelist
mm_miR -> miR
miR[1:100,] -> miR
mm_mRNA -> mRNA
mRNA[1:200,] -> mRNA
AddPrefix(miR, "miR") -> miR_p
AddPrefix(mRNA, "mRNA") -> mRNA_p
#test function
GenesList(method = "s", miR_data = miR_p, mRNA_data = mRNA_p) -> genelist
#check 1
test_that("genelist should be a list of lists", {
expect_equal(class(genelist), "list")
expect_equal(class(genelist[1]), "list")
expect_equal(class(genelist[2]), "list")
expect_equal(class(genelist[1][1]), "list")
})
#internal checks
genedata <- list(miR_data = miR_p, mRNA_data = mRNA_p)
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
expect_equal(genes.split, genelist)
})
#save data
saveRDS(genelist, "genelist_s.rds")

