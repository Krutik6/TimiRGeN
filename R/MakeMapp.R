#' @title MakeMapp
#' @description Create input for the MAPP plugin in pathvisio.
#' @param filt_df Dataframe of mined microRNA-mRNA interactions
#' @param miR_IDs_adj miR_enseml_adj or miR_entrez_adj
#' @param Datatype L for entrez or En for ensembl
#' @return A dataframe which should be saved as a text file for import into
#' pathvisio.
#' @export
#' @usage MakeMapp(filt_df, miR_IDs_adj, Datatype = '')
#' @examples
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' getIDs_miR_mouse(MAE, MAE@ExperimentList$miR) -> MAE
#' 
#' MAE@ExperimentList$Int_matrix <- data.frame(row.names = c(
#' "mmu-miR-320-3p:Acss1","mmu-miR-27a-3p:Odc1"),
#' avecor = c(-0.9191653, 0.7826041),
#' miR = c("mmu-miR-320-3p", "mmu-miR-27a-3p"),
#' mRNA = c("Acss1", "Acss1"),
#' miR_Entrez = c(NA, NA),
#' mRNA_Entrez = c(68738, 18263),
#' TargetScan = c(1, 0),
#' miRDB = c(0, 0),
#' Predicted_Interactions = c(1, 0),
#' miRTarBase = c(0, 1),
#' Pred_Fun = c(1, 1))
#' 
#' MakeMapp(filt_df = MAE@ExperimentList$Int_matrix, 
#' miR_IDs_adj = MAE@ExperimentList$miR_adjusted_entrez, 
#' Datatype = 'L') ->  MAE@ExperimentList$MAPPS
MakeMapp <- function(filt_df, miR_IDs_adj, Datatype){
if (missing(filt_df)) stop('Input miR-mRNA interaction data mined from
MastrixFilter function.');
if (missing(miR_IDs_adj)) stop('Input miR ID data that is adjusted for
repeats. Ensembl or entrez.');
if (missing(Datatype)) stop('En for ensembl or L for entrez.');
merge(x = filt_df, y = miR_IDs_adj, by.x = "miR", by.y = "GENENAME") -> X
if (Datatype == 'L') {
MAPPdata <- data.frame("GENENAME" = X$miR, "ENTREZ" = X$ID,
"code" = 'L', ord = X$mRNA)
} else if(Datatype == 'En'){
MAPPdata <- data.frame("GENENAME" = X$miR, "ENSEMBL" = X$ID,
"code" = 'En', ord = X$mRNA)
} else {print('Enter L or En for data type of interest')}
MAPPdata[order(MAPPdata$ord),] -> MAPPdata
MAPPdata$ord <- NULL
return(MAPPdata)
}
