#' @title linearRegr
#' @description Creates a linear model between a mRNA or miRNA of choice and
#' a selection of it's filtered binding partners. The model is based on
#' the users design.
#' @param mreg Matrix generated from the multiReg function. This should be found as
#' an assay within the MAE used during the multiReg function.
#' @param colselect Column from mreg which contain the gene of interest. Default
#' is 2.
#' @param colpair Column of binding parter of the gene of interest. Default is 3.
#' If NANs are found, test alternative formulas.
#' @param alterpairs Column(s) of other binding partners of the gene of interest.
#' This can include any numer of columns. If NANs are found, test alternative
#' formulas. e.g. = 4:7, c(4, 6, 8), 4.
#' @return A linear regression model which represents miRNA-mRNA interaction(s)
#' which can be further explored.
#' @export
#' @usage linearRegr(mreg, colselect, colpair, alterpairs)
#' @importFrom stats as.formula lm
#' @importFrom stringr str_replace_all
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
#' MAE <- multiReg(MAE = MAE, gene_interest = "Adamts15",
#'                 mRNAreg =TRUE, filt_df=MAE[[11]], miRNA_exp=MAE[[9]],
#'                 mRNA_exp=MAE[[10]])
#'
#' model1 <- linearRegr(mreg = MAE[[12]], colselect =2, colpair =3)
#' summary(model1$regression)
linearRegr <- function(mreg, colselect = 2, colpair = 3, alterpairs){

  if (missing(mreg)) stop('mreg is missing. Add matrix/ assay from the multiReg function.')

  data <- as.data.frame(mreg)

  colnames(data) <- str_replace_all(colnames(data), pattern = "[^[:alnum:]]",
                                    replacement = "")

  vars <- names(data)

  fmla <- paste(vars[colselect], vars[colpair], sep = "~")

  if(!missing(alterpairs)){

    other <- paste(vars[alterpairs], collapse = "+")

    fmla <- paste(fmla, other, sep = "+")

  }

  fmla <- as.formula(fmla)

  list(formula = fmla, regression = lm(fmla, data = data))

}
