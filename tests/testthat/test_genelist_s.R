#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load filtered_genelist
MAE <- readRDS(file = "MAE_Prefix.rds")
miR_p <- assay(MAE, 3)
mRNA_p <- assay(MAE, 4)
#test function
MAE <- genesList(MAE = MAE, method = 's',
                 miR_data = miR_p,
                 mRNA_data = mRNA_p)
#internal checks
genedata <- list(miR_data = miR_p,
                 mRNA_data = mRNA_p)

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
