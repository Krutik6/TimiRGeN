#' @title MakeDynamic
#' @description Make expression data input for pathvisio so changes over time
#' can be visualised.
#' @param miR_expression microRNA data after going through Express function.
#' @param mRNA_expression mRNA data after going through Express function.
#' @param miR_IDs_adj Either miR_ensembl_adj or miR_entrez_adj
#' @param Datatype Either En (ensembl data) or L (entrez data)
#' @return MicroRNA and mRNA dynamic data that can be used in pathvisio if
#'exported.
#' @export
#' @usage MakeDynamic(miR_expression, mRNA_expression, miR_IDs_adj,
#' Datatype = '')
#' @examples
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' library(biomaRt)
#' 
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' 
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' 
#' getIDs_mRNA_mouse(MAE, MAE@ExperimentList$mRNA, "useast") -> MAE
#' getIDs_miR_mouse(MAE, MAE@ExperimentList$miR) -> MAE
#' 
#' Express(df = MAE@ExperimentList$mRNA, dataType = 'Log2FC',
#' genes_ID = MAE@ExperimentList$mRNA_entrez,
#' idColumn = 'GENENAME') -> MAE@ExperimentList$mRNA_express
#' Express(df = MAE@ExperimentList$miR, dataType = 'Log2FC',
#' genes_ID = MAE@ExperimentList$miR_entrez,
#' idColumn = 'GENENAME') -> MAE@ExperimentList$miR_express
#' 
#' MakeDynamic(miR_expression = MAE@ExperimentList$miR_express, 
#' mRNA_expression = MAE@ExperimentList$mRNA_express, 
#' miR_IDs_adj = MAE@ExperimentList$miR_entrez, Datatype = "L"
#' ) -> MAE@ExperimentList$Dynamics
MakeDynamic <- function(miR_expression, mRNA_expression, miR_IDs_adj, Datatype){
if (missing(miR_expression)) stop('Input miR expression data from
smiRk-Express function.');
if (missing(mRNA_expression)) stop('Input mRNA expression data from
smiRk-Express function.');
if (missing(miR_IDs_adj)) stop('Input miR ID data that is adjusted for
repeats. Ensembl or entrez.');
if (missing(Datatype)) stop('En for ensembl or L for entrez.');
rownames(miR_expression) -> miR_expression$names
merge(x= miR_expression, y= miR_IDs_adj, by.x= 'names',
by.y= 'GENENAME', all=TRUE) -> X
rownames(X) <- X$names
X$names <- X$ID.x <- NULL
gsub(colnames(X), pattern = "ID.y", replacement = "ID") -> colnames(X)
names(mRNA_expression) <- names(X)
rbind(X, mRNA_expression) -> Dynamic
cbind(Dynamic, Datatype) -> Dynamic
return(Dynamic)
}
