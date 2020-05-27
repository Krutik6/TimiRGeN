#' @title quickDot
#' @description Creates dotplot which compares the confidence level we have for
#'each wikipathways association to the data, the number of genes associated
#'to each wikipathway and this is all in a ranked fashion for each timepoint
#'(and genetype if s analyses was conducted).
#' @param X Dataframe including count information and wikipathways
#' @param Y Point to the associated name within a list of nested dataframes.
#' @return Dotplot
#' @export
#' @importFrom ggplot2 unit element_rect geom_dotplot
#' @importFrom ggplot2 ggplot scale_fill_continuous labs theme
#' @importFrom ggplot2 element_text ggtitle coord_flip
#' @importFrom stats reorder p.adjust
#' @usage quickDot(X, Y)
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' metadata(MAE)[["e_list"]] <- e_list
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- enrichWiki(MAE = MAE, method = 'c', ID_list = metadata(MAE)[[1]],
#'                    orgDB = org.Mm.eg.db, path_gene = assay(MAE, 3),
#'                    path_name = assay(MAE, 4), ID = "ENTREZID",
#'                    universe = assay(MAE, 3)[[2]])
#'
#' Q <- quickDot(X = metadata(MAE)[[2]][[1]], Y = metadata(MAE)[[2]][1])
quickDot <- function(X, Y){

    if (missing(X)) stop('Input nested dataframe which comes from enrichWiki
                         function. Should be in the metadata as [[]][[]].')

    if (missing(Y)) stop('Input nested dataframe which comes from enrichWiki
                         function. Should be in the metadata as [[]][].')

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
