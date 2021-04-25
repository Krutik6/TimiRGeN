#' @title quickPathwayTC
#' @description Plots time series data from the pathway of interest. Genes
#' which pass a user defined threshold will be highlighted.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param morethan TRUE or FALSE. Default is TRUE.
#' @param threshold Integer that is user defined. Default is 1.
#' @param pathwayname Character which is the name of pathway of interest.
#' Default is "Pathway".
#' @return Line plot of all the genes of interest from a pathway of interest.
#' @export
#' @import gghighlight
#' @usage quickPathwayTC(filt_df, miRNA_exp, mRNA_exp, morethan, threshold,
#'  pathwayname)
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
#' quickPathwayTC(filt_df=MAE[[11]], miRNA_exp=MAE[[9]],
#'                mRNA_exp=MAE[[10]], morethan = TRUE, threshold =1,
#'                pathwayname = "Test")
quickPathwayTC <- function(filt_df, miRNA_exp, mRNA_exp, morethan=TRUE,
                           threshold = 1, pathwayname ="Pathway"){
  if (missing(filt_df)) stop('filt_df is missing. Add assay/ dataframe created by the matrifFilter function.')

  if (missing(miRNA_exp)) stop('miRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on miRNA expression data/ DE data.')

  if (missing(mRNA_exp)) stop('mRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on mRNA expression data/ DE data.')

  Time <- Expression <- Gene <-

  Prep <- clustPrep(filt_df, miRNA_exp, mRNA_exp)

  if (morethan==TRUE) {

    ggplot(Prep, aes(Time, Expression, color=Gene)) +

      geom_line(stat="identity", size = 1.5)+

      theme_bw() +

      theme(legend.position = 'none')+

      labs(title= paste0("Genes across ", pathwayname),

           x="Time",

           y="Scaled Expression")+

      theme(plot.title=element_text(size=20, face="bold",hjust = 0.5),

            axis.text.x=element_text(size=15),

            axis.text.y=element_text(size=15),

            axis.title.x=element_text(size=17),

            axis.title.y=element_text(size=17))+

      theme(axis.line = element_line(colour = "black"),

            panel.grid.major = element_blank(),

            panel.grid.minor = element_blank(),

            panel.border = element_blank(),

            panel.background = element_blank())+

      gghighlight(max(Expression) > threshold,

                  max_highlight = 1,

                  use_direct_label = TRUE)

  }else if (morethan==FALSE) {

    ggplot(Prep, aes(Time, Expression, color=Gene)) +

      geom_line(stat="identity", size = 1.5)+

      theme_bw() +

      theme(legend.position = 'none')+

      labs(title=paste0("Genes across ", pathwayname),

           x="Time",

           y="Scaled Expression")+

      theme(plot.title=element_text(size=20, face="bold",hjust = 0.5),

            axis.text.x=element_text(size=15),

            axis.text.y=element_text(size=15),

            axis.title.x=element_text(size=17),

            axis.title.y=element_text(size=17))+

      theme(axis.line = element_line(colour = "black"),

            panel.grid.major = element_blank(),

            panel.grid.minor = element_blank(),

            panel.border = element_blank(),

            panel.background = element_blank())+

      gghighlight(max(Expression) < threshold, max_highlight = 1,

                  use_direct_label = TRUE)

  }else{print("morethan parameter must be TRUE or FALSE")

  }
}
