#' @title quickDMap
#' @description Creates a companion heatmap for the dendrogram made by
#' quickDendro.
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
#' @return A heatmap with time points as the x axis and genes as the y axis.
#' Gene order will be the same as quickDendro.
#' @export
#' @usage quickDMap(filt_df, miRNA_exp, mRNA_exp, distmeth, hclustmeth,
#' pathwayname)
#' @examples
#' library(org.Mm.eg.db)
#' miR <- mm_miR[1:100,]
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus')
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
#'
#' quickDMap(filt_df=MAE[[11]], miRNA_exp=MAE[[9]],
#'           mRNA_exp=MAE[[10]], pathwayname = "Test")
quickDMap <- function(filt_df, miRNA_exp, mRNA_exp, distmeth="maximum",
                      hclustmeth = "ward.D", pathwayname = "Pathway"){

  Gene <- Expression <- tracecol <- NULL

  fit <- hClustPrep(filt_df, miRNA_exp, mRNA_exp, distmeth,
                    hclustmeth)

  D <- ggdendrogram(fit, rotate = TRUE, theme_dendro = FALSE) +
    theme_bw() +
    labs(title= paste("Gene clusters within", pathwayname),
         x="Genes",
         y="Distance")+
    theme(plot.title=element_text(size=15, face="bold",hjust = -8),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=17),
          axis.title.y=element_text(size=17))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

  labs <- D$plot_env$data$labels

  olabs <- labs[order(labs$x, decreasing = TRUE),]

  Prep <- clustPrep(filt_df, miRNA_exp, mRNA_exp)

  X <- Prep %>% spread(Gene, Expression)

  rownames(X) <- X$Time

  X$Time <- NULL

  Y <- as.matrix(t(X))

  Z <- Y[match(olabs$label, rownames(Y)), ]

  my_palette <- grDevices::colorRampPalette(c("white", "cyan", "violet"))(n = 200)

  par(oma =c (1,1,0,5))

  heatmap.2(Z, trace = "n", col = my_palette,Colv = FALSE, dendrogram = "none",
            Rowv=FALSE, key.xlab = "Expression", key.title = "", key.ylab = "",
            ylab = "", xlab = "Time",
            main = paste0(pathwayname ," Gene Expression"), cexRow = 1.1,
            cexCol = 2,  keysize=.8, key.par = list(cex=.7),
            density.info=c("none"), denscol=tracecol)
}
