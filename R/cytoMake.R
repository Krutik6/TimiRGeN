#' @title cytoMake
#' @description Creates a cytoscape network based on the output of matrixFilter.
#' Requires cytoscapePing() to be used. Make sure Cytoscape is open first. Must
#' use Cytoscape version 3.7 or later.
#' @param interactionData Dataframe which contains filtered miR-mRNA
#' interactions. This is output from matrixFilter and should be found as an
#' assay within the MAE used in the matrixFilter function.
#' @param titleString Title of the network. Enter a string which cytoscape
#' will see as the graph name. Default is "Network".
#' @param collectionString Title of the collection of networks. Enter string
#' which cytoscape will see as the collection name. Many differently titled
#' graphs can be added to a single collection. Default is
#' "miR-mRNA interactions".
#' @return A network of filtered miR-mRNA interactions specific for a pathway
#' of interest. It will be visible in cytoscape version 3.7 or later.
#' @export
#' @importFrom RCy3 createNetworkFromDataFrames setVisualStyle layoutNetwork
#' @importFrom RCy3 cytoscapePing
#' @usage cytoMake(interactionData, titleString = '', collectionString = '')
#' @examples
#' \dontrun{
#' Filt_df <- data.frame(row.names = c("mmu-miR-320-3p:Acss1",
#'                                     "mmu-miR-27a-3p:Odc1"),
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
#' RCy3::cytoscapePing()
#'
#' cytoMake(interactionData = Filt_df, titleString = 'test' ,
#'          collectionString = 'collectiontest')
#' }
cytoMake <- function(interactionData, titleString = "Network",
                     collectionString = "miR-mRNA interactions"){

    if (missing(interactionData)) stop('interactionData is missing. Add filtered miR-mRNA dataframe. Please use the matrixFilter function first. Output of the matrixFilter function should be stored as an assay within the MAE used in the matrixFilter function.')

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
    RCy3::createNetworkFromDataFrames(nodes,edges,
                                title= titleString,
                                collection= collectionString)

    # Edit visuals
    RCy3::setVisualStyle('Marquee')

    RCy3::layoutNetwork('cose defaultSpringCoefficient=0.000008
                         defaultSpringLength=50')
}
