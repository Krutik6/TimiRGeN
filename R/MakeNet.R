#' @title MakeNet
#' @description Creates an igraph object from filtered miR-mRNA interactions
#' from database mining.
#' @param filt_df Filtered miR-mRNA interactions from the mining matrix
#' produced by the DataMiningMatrix function
#' @return A list which is the input for Quicknet
#' @export
#' @importFrom igraph graph_from_data_frame
#' @usage MakeNet(filt_df)
#' @examples
#' filt_df <- structure(list(avecor = c(-0.929199786400515, -0.729228501795928,
#'-0.431983639087243, -0.55088842103792, -0.978422379116014,
#'-0.627856061946295, -0.998864281242427, -0.877362389481312,
#'-0.990125035397228,
#'-0.771338310408749), miR = structure(c(9L, 5L, 6L, 2L, 8L, 4L, 3L, 4L,
#'7L, 1L), .Label = c("hsa-miR-107", "hsa-miR-193a-3p", "hsa-miR-28-5p",
#'"hsa-miR-331-3p", "hsa-miR-362-3p", "hsa-miR-362-5p", "hsa-miR-429",
#'"hsa-miR-590-5p", "hsa-miR-630" ), class = "factor"),
#' mRNA = structure(c(1L, 2L, 2L, 3L, 3L, 4L, 5L, 5L,
#' 5L, 6L), .Label = c("IGF1R", "PRKCA", "TESK2", "THBS1",  "TLN2", "VAV3"),
#'class = "factor"), miR_Entrez = structure(c(8L,  7L, 7L, 4L, 6L, 3L, 5L,
#' 3L, 1L, 2L), .Label = c("ENSG00000198976",  "ENSG00000198997",
#'"ENSG00000199172", "ENSG00000207614", "ENSG00000207651",
#'"ENSG00000207741", "ENSG00000208015", "ENSG00000283798"),
#'class = "factor"), mRNA_Entrez = structure(c(4L, 5L, 5L, 1L, 1L, 3L, 6L,
#'6L, 6L, 2L), .Label = c("ENSG00000070759", "ENSG00000134215",
#'"ENSG00000137801", "ENSG00000140443", "ENSG00000154229","ENSG00000171914"),
#'class = "factor"), TargetScan = c(0L,  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
#' 1L), miRDB = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
#' Predicted_Interactions = c(1L,  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
#' miRTarBase = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), Pred_Fun = c(2L,
#' 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)), class = "data.frame",
#' row.names = c("hsa-miR-630:IGF1R", "hsa-miR-362-3p:PRKCA",
#'"hsa-miR-362-5p:PRKCA", "hsa-miR-193a-3p:TESK2", "hsa-miR-590-5p:TESK2",
#'"hsa-miR-331-3p:THBS1", "hsa-miR-28-5p:TLN2", "hsa-miR-331-3p:TLN2",
#'"hsa-miR-429:TLN2", "hsa-miR-107:VAV3" ))
#' MakeNet(filt_df = filt_df) -> net
MakeNet <- function(filt_df){
rbind(filt_df, filt_df) -> df
df$id <- "id"
genes <- data.frame(genes = c(as.character(filt_df$miR),
as.character(filt_df$mRNA)))
as.integer(factor(genes[,1])) -> genes$id
paste("s", genes$id, sep = "") -> df$id
rownames(df)<- NULL
max(as.integer(rownames(df))/2) -> halfway
nodes <- data.frame(id = df$id,
genes = c(as.character(df$miR[seq_len(halfway)]),
as.character(df$mRNA[seq_len(halfway)])))
nodes$genetype <- "mRNA"
max(as.integer(rownames(nodes))/2) -> number_miRs
for (i in seq_len(number_miRs)) {
nodes$genetype[i] <- "miR"
}
links <- data.frame(from = nodes$id[seq_len(halfway)],
to = nodes$id[-c(seq_len(halfway))],
Databases = filt_df$Pred_Fun,
Correlation = filt_df$avecor,
type = "hyperlink")
nodes[! duplicated(nodes$genes),] -> nodes
net <- graph_from_data_frame(d=links, vertices=nodes, directed=TRUE)
return(net)
}
