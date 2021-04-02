#' @title quickTC
#' @description Plots miRNA:mRNA pair over timecourse.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param pair Interger representing the pair to be explored.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data
#' @param scale TRUE or FALSE. Should data be scales. Default is FALSE.
#' @param Interpolation TRUE or FALSE. Should the whole time course be
#' interpolated over by a smooth spline? Default is FALSE. This is most useful
#' for longer time courses.
#' @param timecourse If Iterpolation is TRUE, how many time points should be
#' interpolated over?
#' @return Time course plot of selected pair.
#' @export
#' @usage quickTC(filt_df, pair, miRNA_exp, mRNA_exp, scale,Interpolation,
#' timecourse)
#' @importFrom ggplot2 theme_classic scale_colour_manual
#' @importFrom FreqProf approxm
#' @examples
#' library(org.Mm.eg.db)
#' 
#' miR <- mm_miR[1:100,]
#' 
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#' 
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
#' quickTC(filt_df=MAE[[11]], pair=1, miRNA_exp=MAE[[9]],
#'         mRNA_exp=MAE[[10]], scale = FALSE)
quickTC <- function(filt_df, pair, miRNA_exp, mRNA_exp, scale=FALSE,
                    Interpolation=FALSE, timecourse){

  Time <- Expression <- Gene <- NULL

  Int <- pickPair(filt_df, pair, miRNA_exp, mRNA_exp, scale)

  if (Interpolation == TRUE) {

    if (missing(timecourse)) stop('timecourse is missing. How many time points to interpolate over? This should be the whole time course.')

    Int <- FreqProf::approxm(as.data.frame(Int), timecourse,

                             method = "spline")

  }else if (Interpolation == FALSE) {

    Int <- Int

  }

  Ranks <- filt_df[,c(1,2,3)][order(filt_df$corr, decreasing = FALSE),]

  Corr <- round(Ranks[1][pair,],2)

  X <- as.data.frame(Int)

  rownames(X) <- X$Time

  X$Time <- NULL

  Melted <- melt(as.matrix(X), varnames = c("Time", "Gene"))

  colnames(Melted)[3] <- "Expression"

  miR <- unique(Melted$Gene)[[1]]

  mRNA <- unique(Melted$Gene)[[2]]

  ggplot(Melted, aes(x = Time, y = Expression, group = Gene, color = Gene)) +

    geom_line(data = ~ subset(., Gene == paste(miR)), size =2) +

    geom_line(data = ~ subset(., Gene == paste(mRNA)), size = 2) +

    scale_colour_manual(values=c("Red", "Blue"))+

    theme_classic()+

    labs(title= paste0(miR, ":", mRNA, " Expression"),

         x="Time",

         y="Expression",

         subtitle=paste0("Corr = ", Corr))+

    theme(plot.title=element_text(size=20, face="bold",hjust = 0.5),

          axis.text.x=element_text(size=15),

          axis.text.y=element_text(size=15),

          axis.title.x=element_text(size=20),

          axis.title.y=element_text(size=20),

          legend.text=element_text(size=12))+

    theme(plot.subtitle=element_text(size=25, hjust=1.2,

                                     face="italic", color="black"))
}
