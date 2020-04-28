#' @title cytoMake
#' @description Creates a cytoscape network based on the output of MatrixFilter.
#'Requires RCy3 package and cytoscapePing() to be used. Make sure Cytoscape
#'is open first.
#' @param interactionData Output from MatrixFiler. Requires a dataframe with
#'columns names MIR and mRNA
#' @importFrom RCy3 createNetworkFromDataFrames setVisualStyle layoutNetwork
#' @param titleString Title of the network.
#' @param collectionString Title of the collection fo networks.
#' @return A network visible in cytoscape.
#' @export
#' @usage cytoMake(interactionData, titleString = '', collectionString = '')
#' @examples
#' \donttest{
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#' MAE <- getIDs_miR_mouse(MAE, assay(MAE, 1))
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
#' cytoscapePing()
#'
#' cytoMake(interactionData = Filt_df, titleString = 'test' ,
#'          collectionString = 'collectiontest')
#' }
cytoMake <- function(interactionData, titleString, collectionString){

    if (missing(interactionData)) stop('Add filtered miR-mRNA dataframe');
    if (missing(titleString)) stop('Add title of network');
    if (missing(collectionString)) stop('Add title of network collections');

    interaction_data <- interactionData
    # Create nodes
    nodes <- data.frame(id = c(as.character(interaction_data$miR),
                               as.character(interaction_data$mRNA)),
                        stringsAsFactors = FALSE)
    # Create edges
    edges <- data.frame(source = interaction_data$miR,
                        target = interaction_data$mRNA,
                        interaction = 'inhibits',
                        weight = NA,
                        stringsAsFactors = FALSE)
    # Create network
    createNetworkFromDataFrames(nodes,edges,
                                title= titleString,
                                collection= collectionString)
    # Edit visuals
    setVisualStyle('Marquee')

    layoutNetwork('cose defaultSpringCoefficient=0.000008
                  defaultSpringLength=50')
}
