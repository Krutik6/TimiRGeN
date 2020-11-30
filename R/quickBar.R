#' @title quickBar
#' @description Creates a bar plot which compares the confidence levels for
#' each wikipathway association to the filtered input data. The number of genes
#' in common between a pathway and the input data are taken into account to
#' generate confidence scores. This function is used specifically for a single
#' time point at a time, or if using "s" analysis, the function is used for a
#' single time point within a single gene type (miR or mRNA).
#' @param X Dataframes within a list including count information, confidence
#' scores and wikipathway information. This is the output from the enrichWiki
#' function. It will be stored as metadata within the MAE used in the enrichWiki
#' function. Data can be retrieved using [[i]][[j]] on the output of enrichWiki.
#' @param Y String which is associated to the nested dataframe selected for X.
#' This is the output from the enrichWiki function. It will be stored as
#' metadata within the MAE used in the enrichWiki function. Data can
#' be retrieved using [[i]][j] on the output of enrichWiki.
#' @return Bar plot showing which pathways are most enriched for genes found
#' at each time point ("c") or at each time point within a gentype ("s").
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_continuous labs theme
#' @importFrom  ggplot2 element_text ggtitle coord_flip
#' @importFrom stats reorder p.adjust
#' @usage quickBar(X, Y)
#' @examples
#' library(org.Mm.eg.db)
#'
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list_mouse
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- enrichWiki(MAE = MAE, method = 'c', ID_list = metadata(MAE)[[1]],
#'                    orgDB = org.Mm.eg.db, path_gene = assay(MAE, 1),
#'                    path_name = assay(MAE, 2), ID = "ENTREZID",
#'                    universe = assay(MAE, 1)[[2]])
#'
#' q <- quickBar(X = metadata(MAE)[[2]][[1]], Y = metadata(MAE)[[2]][1])
#'
#' # to view bar plot enter plot(q)
quickBar <- function(X, Y){

    if (missing(X)) stop('
                         X is missing.
                         Input nested dataframe which is output from
                         enrichWiki. Should be in the metadata of the MAE used
                         in the enrichWiki function, and this can be retrieved
                         using metadata(MAE)[[i]][[j]].')

    if (missing(Y)) stop('
                         Y is missing.
                         Input name of the nested dataframe associated with X.
                         This is output from enrichWiki. Should be in the
                         metadata of the MAE used in the enrichWiki function,
                         and this can be retrieved using metadata(MAE)[[i]][j].
                         ')

  Description <- Count <- NULL

    ggplot2::ggplot(head(X[which(X$p.adjust < 0.05),],n = 15),
                    aes(x=reorder(Description, -p.adjust),
                        y=Count,
                        fill=-p.adjust)) +

      geom_bar(stat = "identity",
               width = 0.5) +

      scale_fill_continuous(type = "gradient") +

      labs(y = "Associated genes",
           x = "Wikipathways",
           fill = "p.adjust") +

      ggtitle(names(Y)) +

      theme(plot.title = element_text(size=30, face = "bold", hjust = 0.5),
            legend.key.size = unit(5, "line"),
            legend.text=element_text(size=15),
            legend.title =element_text(size=20),
            axis.text=element_text(size=30),
            axis.title = element_text(size = 20)) +

      theme(panel.background = element_rect(fill = 'white',
                                            colour = 'black'))+
      coord_flip()
}
