#' @title quickDot
#' @description Creates a dot plot which compares the confidence levels for
#' each wikipathway association to the data. The number of genes in common
#' between a pathway and the input data are taken into account to generate
#' confidence scores. This function is to be used specifically for a single
#' time point at a time, or if s analysis has been done, then a single time
#' point within a single gene type (mir or mRNA).
#' @param X Dataframe within a list including count information,confidence
#' scores and wikipathway information. This is the output from the enrichWiki
#' function. It will be stored as metadata within the MAE used in the enrichWiki
#' function. Data can be retrieved using [[i]][[j]] on the output of enrichWiki.
#' @param Y String which is associated to the name of a nested
#' dataframes. This is the output from the enrichWiki function. It will be
#' stored as metadata within the MAE used in the enrichWiki function. Data can
#' be retrieved using [[i]][j] on the output of enrichWiki.
#' @return Dot plot showing which pathways are most enriched for genes found
#' in a time point found in all gene types or in a time point within a gene
#' type.
#' @export
#' @importFrom ggplot2 unit element_rect geom_dotplot
#' @importFrom ggplot2 ggplot scale_fill_continuous labs theme
#' @importFrom ggplot2 element_text ggtitle coord_flip
#' @importFrom stats reorder p.adjust
#' @usage quickDot(X, Y)
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
#' q <- quickDot(X = metadata(MAE)[[2]][[1]], Y = metadata(MAE)[[2]][1])
#'
#' # to view dot plot enter plot(q)
quickDot <- function(X, Y){

    if (missing(X)) stop('Input nested dataframe which is output from
                           enrichWiki. Should be in the metadata of the MAE used
                           in the enrichWiki function, and can
                           be retrieved using [[i]][[j]].')

    if (missing(Y)) stop('Input name of the  nested dataframe which is output
                           from enrichWiki. Should be in the metadata of the MAE
                           used in the enrichWiki function, and can be retrieved
                           using [[i]][j].')

    Description <- Count <- NULL

    ggplot2::ggplot()+

    geom_dotplot(mapping = aes(x= reorder(Description, -p.adjust),
                               y=Count,
                               fill=-p.adjust),
                 data = head(X[which(X$p.adjust < 0.05),], n = 15),
                 binaxis = 'y',
                 dotsize = 2,
                 method = 'dotdensity',
                 binpositions = 'bygroup',
                 binwidth = 0.2,
                 stackdir = "center") +

    scale_fill_continuous(type = "gradient") +

    labs(y = "Associated genes",
         x = "wikipathways",
         fill = "p.adjust") +

    theme(axis.text=element_text(size=16)) +

    ggtitle(names(Y)) +

    theme(plot.title = element_text(2, face = "bold", hjust = 0.5),
          legend.key.size = unit(2, "line")) +

    theme(panel.background = element_rect(fill = 'white', colour = 'black'))+

    coord_flip()
    }
