#' @title MakeDynamic
#' @description Make expression data input for pathvisio so changes over time
#' can be visualised.
#' @param MAE MultiAssayObject.
#' @param miR_expression microRNA data after going through Express function.
#' @param mRNA_expression mRNA data after going through Express function.
#' @param miR_IDs_adj Either miR_ensembl_adj or miR_entrez_adj
#' @param dataType Either En (ensembl data) or L (entrez data)
#' @return MicroRNA and mRNA dynamic data that can be used in pathvisio if
#'exported.
#' @export
#' @usage MakeDynamic(MAE, miR_expression, mRNA_expression, miR_IDs_adj,
#'                    dataType = '')
#' @examples
#' library(biomaRt)
#' library(MultiAssayExperiment)
#' miR <- mm_miR[1:100,]
#' mRNA <- mm_mRNA[1:200,]
#' MAE <- StartObject(miR = miR, mRNA = mRNA)
#' MAE <- getIDs_miR_mouse(MAE, assay(MAE, 1))
#' MAE <- getIDs_mRNA_mouse(MAE, assay(MAE, 2), "useast")
#' MAE <- Express(MAE, df = assay(MAE, 1), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 3),
#'                idColumn = 'GENENAME',
#'                name = "miR_express")
#' MAE <- Express(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 7),
#'                idColumn = 'GENENAME',
#'                name = 'mRNA_express')
#' 
#' MAE <- MakeDynamic(MAE, miR_expression = assay(MAE, 9), 
#'                   mRNA_expression = assay(MAE, 10), 
#'                   miR_IDs_adj = assay(MAE, 5), 
#'                   dataType = "L")
MakeDynamic <- function(MAE, miR_expression, mRNA_expression, miR_IDs_adj, 
                        dataType){
    
    if (missing(MAE)) stop('Add MAE object');
    if (missing(miR_expression)) stop('Input miR expression data from
    smiRk-Express function.');
    if (missing(mRNA_expression)) stop('Input mRNA expression data from
    smiRk-Express function.');
    if (missing(miR_IDs_adj)) stop('Input miR ID data that is adjusted for
    repeats. Ensembl or entrez.');
    if (missing(dataType)) stop('En for ensembl or L for entrez.');
    
    miR_expression <- miR_expression
    mRNA_expression <- mRNA_expression
    miR_IDs_adj <- miR_IDs_adj
    
    miR_expression$names <- rownames(miR_expression)
    # Merge time series data with IDs
    X <- merge(x= miR_expression, y= miR_IDs_adj, by.x= 'names',
               by.y= 'GENENAME', all=TRUE)
    
    rownames(X) <- X$names
    X$names <- X$ID.x <- NULL
    
    colnames(X) <- gsub(colnames(X), pattern = "ID.y", replacement = "ID") 
    
    names(mRNA_expression) <- names(X)
    
    Dynamic <- rbind(X, mRNA_expression)
    Dynamic <- cbind(Dynamic, dataType)
    
    MAE2 <- suppressMessages(MultiAssayExperiment(list('Dynamic' = Dynamic)))
    MAE <- c(MAE, MAE2)
    
return(MAE)
}
