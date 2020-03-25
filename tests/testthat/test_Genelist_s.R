#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(stringr)
#load filtered_genelist
MAE <- readRDS(file = "MAE_Prefix.rds")

#test function
MAE <- GenesList(MAE = MAE, method = 's', 
                 miR_data = assay(MAE, 3), 
                 mRNA_data = assay(MAE, 4))

#check 1
test_that("genelist should be a list of lists", {
    expect_equal(class(metadata(MAE)[1]), "list")
})
#internal checks
genedata <- list(miR_data = assay(MAE, 3),
                 mRNA_data = assay(MAE, 4))
#check 2
test_that("genedata is a list", {
    expect_equal(class(genedata), "list")
})
#contineu
geneslist <- lapply(genedata, function(df)
    {
    Unigenes <- unique(str_extract(names(df),"\\S+\\."))
    List <- lapply(Unigenes,function(name){return(df[,grep(name,names(df),
                                                           fixed=TRUE)])})
    names(List) <- Unigenes
return(List)
})
metadata(MAE)[["geneslist"]] <- geneslist
metadata(MAE)
#save data
saveRDS(metadata(MAE)[[1]], "genelist_s.rds")