#' @title CytoMake
#' @description Creates a cytoscape network based on the output of MatrixFilter.
#'Requires RCy3 package and cytoscapePing() to be used. Make sure Cytoscape
#'is open first.
#' @param interaction_data Output from MatrixFiler. Requires a dataframe with
#'columns names MIR and mRNA
#' @importFrom RCy3 createNetworkFromDataFrames setVisualStyle layoutNetwork
#' @param titleString Title of the network.
#' @param collectionString Title of the collection fo networks.
#' @return A network visible in cytoscape.
#' @export
#' @usage CytoMake(interaction_data, titleString = '', collectionString = '')
CytoMake <- function(interaction_data, titleString, collectionString){
nodes <- data.frame(id = c(as.character(interaction_data$miR),
as.character(interaction_data$mRNA)),
stringsAsFactors = FALSE)
edges <- data.frame(source = interaction_data$miR,
target = interaction_data$mRNA,
interaction = 'inhibits',
weight = NA,
stringsAsFactors = FALSE)
createNetworkFromDataFrames(nodes,edges, title= titleString,
collection= collectionString)
setVisualStyle('Marquee')
layoutNetwork('cose defaultSpringCoefficient=0.000008
defaultSpringLength=50')
}
