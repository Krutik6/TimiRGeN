#' @title makeMapp
#' @description Creates input for the MAPP plugin in pathvisio. This will
#' help to add the filtered miRNAs to the wikipathway of interest.
#' Follow vignette instructions on how to save this file and further
#' instructions found in /inst/Pathvisio_GRN_guide.pdf to see how
#' this can help in GRN construction.
#' @param MAE MultiAssayExperiment object to store the output of makeMapp. It
#' is recommended to use the same MAE object which stores the results from
#' matrixFilter.
#' @param filt_df Dataframe of mined microRNA-mRNA interactions. This is output
#' of matrixFilter. It should be stored as an assay in the MAE used in the
#' matrixFilter function.
#' @param miR_IDs_adj Dataframes with adjusted gene IDs to account for -5p and
#' -3p specific miRs. miR_adjusted_entrez or miR_adjusted_ensembl. Should be
#' found as an assay in the MAE used a getIdsMir function.
#' @param dataType String which represents the gene ID used in this analysis.
#' Either "En" (ensembl data) or "L" (entrez data).
#' @return A dataframe containing microRNA and adjusted gene ID information
#' which should be saved as a text file for import into pathvisio via the MAPP
#' app.
#' @export
#' @usage makeMapp(MAE, filt_df, miR_IDs_adj, dataType = '')
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
#'
#' Filt_df <- data.frame(row.names = c("mmu-miR-320-3p:Acss1",
#'                                      "mmu-miR-27a-3p:Odc1"),
#'                       avecor = c(-0.9191653, 0.7826041),
#'                       miR = c("mmu-miR-320-3p", "mmu-miR-27a-3p"),
#'                       mRNA = c("Acss1", "Acss1"),
#'                       miR_Entrez = c(NA, NA),
#'                       mRNA_Entrez = c(68738, 18263),
#'                       TargetScan = c(1, 0),
#'                       miRDB = c(0, 0),
#'                       Predicted_Interactions = c(1, 0),
#'                       miRTarBase = c(0, 1),
#'                       Pred_Fun = c(1, 1))
#'
#' MAE <- makeMapp(MAE, filt_df = Filt_df, miR_IDs_adj = assay(MAE, 5),
#'                 dataType = 'L')
makeMapp <- function(MAE, filt_df, miR_IDs_adj, dataType){

    if (missing(MAE)) stop('Add MAE object. This will store the output of
                           makeMapp Please use matrixFilter first.')

    if (missing(filt_df)) stop('Add dataframe of filtered miR-mRNA
                               interactions. Please use the matrixFilter
                               function first. Output of matrixFilter should
                               be stored as an assay within the MAE used in the
                               matrixFilter function.')

    if (missing(miR_IDs_adj)) stop('Add adjusted miR gene ID data.
                                   Please use the getIdsMirHuman or
                                   getIdsMirMouse function first. Output of
                                   a getIdsMir function should be stored as an
                                   assay within the MAE used in the
                                   getIdsMir function.')

    if (missing(dataType)) stop('Add a string. "En" for ensembl or "L"
                                for entrez.')

    # Merge filtered data frame of interactions and adjusted miR IDs
    X <- merge(x = filt_df, y = miR_IDs_adj, by.x = "miR", by.y = "GENENAME")

    # For entrezgene IDs
    if (dataType == 'L') {

        MAPPdata <- data.frame("GENENAME" = X$miR, "ENTREZ" = X$ID,
                               "code" = 'L', ord = X$mRNA)
    # For ensmbl gene IDs
    } else if(dataType == 'En'){

        MAPPdata <- data.frame("GENENAME" = X$miR, "ENSEMBL" = X$ID,
                               "code" = 'En', ord = X$mRNA)

    } else {print('Enter L or En for data type of interest')}

    MAPPdata <- MAPPdata[order(MAPPdata$ord),]

    MAPPdata$ord <- NULL

    MAE2 <- suppressMessages(MultiAssayExperiment(list('MAPPdata' = MAPPdata)))

    MAE <- c(MAE, MAE2)

return(MAE)
}
