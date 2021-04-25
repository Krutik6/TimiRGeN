#' @title quickDendro
#' @description Ceates a dendrogram of the genes from the pathway of interest.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param distmeth Dist method for hierarchical clustering. Default is
#' "maximum".
#' @param hclustmeth Hclust method for hierarchical clustering. Default is
#' "ward.D".
#' @param pathwayname Character which is the name of pathway of interest.
#' Default is "Pathway".
#' @return A dendrogram with the genes on the Y axis and the distances on the
#' X axis.
#' @export
#' @usage quickDendro(filt_df, miRNA_exp, mRNA_exp, distmeth, hclustmeth,
#' pathwayname)
#' @importFrom ggdendro ggdendrogram
#' @importFrom ggplot2 theme_bw element_line element_blank
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR[1:50,]
#'
#' mRNA <- mm_mRNA[1:100,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus', orgDB = org.Mm.eg.db)
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 1), dataType = 'Log2FC',
#'                      genes_ID = assay(MAE, 3),
#'                       idColumn = 'GENENAME',
#'                       name = "miRNA_log2fc")
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                      genes_ID = assay(MAE, 7),
#'                      idColumn = 'GENENAME',
#'                      name = "mRNA_log2fc")
#'
#' Filt_df <- data.frame(row.names = c("mmu-miR-145a-3p:Adamts15",
#'                                    "mmu-miR-146a-5p:Acy1"),
#'                      corr = c(-0.9191653, 0.7826041),
#'                      miR = c("mmu-miR-145a-3p", "mmu-miR-146a-5p"),
#'                      mRNA = c("Adamts15", "Acy1"),
#'                      miR_Entrez = c(387163, NA),
#'                      mRNA_Entrez = c(235130, 109652),
#'                      TargetScan = c(1, 0),
#'                      miRDB = c(0, 0),
#'                      Predicted_Interactions = c(1, 0),
#'                      miRTarBase = c(0, 1),
#'                      Pred_Fun = c(1, 1))
#'
#' MAE <- matrixFilter(MAE, miningMatrix = Filt_df, negativeOnly = FALSE,
#'                    threshold = 1, predictedOnly = FALSE)
#'
#' quickDendro(filt_df=MAE[[11]], miRNA_exp=MAE[[9]],
#'             mRNA_exp=MAE[[10]], pathwayname = "Test")
quickDendro <- function(filt_df, miRNA_exp, mRNA_exp, distmeth="maximum",
                       hclustmeth = "ward.D", pathwayname = "Pathway"){

  fit <- hClustPrep(filt_df, miRNA_exp, mRNA_exp, distmeth, hclustmeth)

  ggdendrogram(fit, rotate = TRUE, theme_dendro = FALSE) +

    ggplot2::theme_bw() +

    labs(title= paste("Gene clusters within", pathwayname),

         x="Genes",

         y="Distance")+

    theme(plot.title=element_text(size=15, face="bold",hjust = 0),

          axis.text.x=element_text(size=15),

          axis.text.y=element_text(size=12),

          axis.title.x=element_text(size=17),

          axis.title.y=element_text(size=17))+

    theme(axis.line = element_line(colour = "black"),

          panel.grid.major = element_blank(),

          panel.grid.minor = element_blank(),

          panel.border = element_blank(),

          panel.background = element_blank())

}
