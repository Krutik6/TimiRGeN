#' @title Quickdot
#' @description Creates dotplot which compares the confidence level we have for
#'each wikipathways association to the data, the number of genes associated
#'to each wikipathway and this is all in a ranked fashion for each timepoint
#'and genetype.
#' @param X Dataframe including count information and wikipathways
#' @param Y Point to the associated name within a list of nested dataframes.
#' @return dotplot
#' @export
#' @importFrom ggplot2 unit element_rect geom_dotplot
#' @usage Quickdot(X, Y)
#' @examples
#' library(ggplot2)
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' e_list -> elist
#' downloadGMT(speciesInitial = "Mm")
#' EnrichWiki(method = 'c', elist, org.Mm.eg.db,
#' path_gene, path_name, ID = 'ENTREZID',
#' universe = path_gene$gene) -> sigwiki
#' Quickdot(X = sigwiki$D7_wikipathways, Y = sigwiki[4]) -> Q
Quickdot <- function(X, Y){
if (missing(X)) stop('Input nested dataframe from list.');
if (missing(Y)) stop('Input name of nested dataframe');
Description <- Count <- NULL
ggplot()+
geom_dotplot(mapping = aes(x= reorder(Description, -p.adjust), y=Count,
fill=-p.adjust),
data = head(X[which(X$p.adjust < 0.05),], n = 15),
binaxis = 'y', dotsize = 2,
method = 'dotdensity', binpositions = 'bygroup', binwidth = 0.3,
stackdir = "center")+
scale_fill_continuous(type = "gradient") +
labs(y = "Associated genes", x = "wikipathways", fill = "p.adjust") +
theme(axis.text=element_text(size=14)) +
ggtitle(names(Y)) +
theme(plot.title = element_text(2, face = "bold", hjust = 0.5),
legend.key.size = unit(2, "line")) +
theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
coord_flip()
}
