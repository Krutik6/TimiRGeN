#' @title makeNet
#' @description Creates an igraph object from filtered miR-mRNA interactions
#' from database mining.
#' @param MAE MultiAssayExperiment Object.
#' @param filt_df Filtered miR-mRNA interactions from the mining matrix
#' produced by the DataMiningMatrix function.
#' @return A list which is the input for Quicknet
#' @export
#' @importFrom igraph graph_from_data_frame
#' @usage makeNet(MAE, filt_df)
#' @examples
#' library(igraph)
#' library(biomaRt)
#' library(MultiAssayExperiment)
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
#'MAE <- MultiAssayExperiment()
#'MAE <- makeNet(MAE, Filt_df)
makeNet <- function(MAE, filt_df){

    if (missing(MAE)) stop('Add MAE object');
    if (missing(filt_df)) stop('Add filtered miR-mRNA interaction data frame')

    metadata <- `metadata<-` <- NULL

    df <- rbind(filt_df, filt_df)
    df$id <- "id"
    genes <- data.frame(genes = c(as.character(filt_df$miR),
                        as.character(filt_df$mRNA)))
    genes$id<- as.integer(factor(genes[,1]))
    df$id <- paste("s", genes$id, sep = "")
    rownames(df) <- NULL
    halfway <- max(as.integer(rownames(df))/2)

    nodes <- data.frame(id = df$id,
                        genes = c(as.character(df$miR[seq_len(halfway)]),
                        as.character(df$mRNA[seq_len(halfway)])))
    nodes$genetype <- "mRNA"
    number_miRs <- max(as.integer(rownames(nodes))/2)
    for (i in seq_len(number_miRs)) {
        nodes$genetype[i] <- "miR"
    }

    links <- data.frame(from = nodes$id[seq_len(halfway)],
                        to = nodes$id[-c(seq_len(halfway))],
                        Databases = filt_df$Pred_Fun,
                        Correlation = filt_df$avecor,
                        type = "hyperlink")

    nodes <- nodes[! duplicated(nodes$genes),]
    net <- graph_from_data_frame(d=links, vertices=nodes, directed=TRUE)
    metadata(MAE)[["net"]] <- net
return(MAE)
}
