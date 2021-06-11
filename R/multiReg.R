#' @title multiReg
#' @description Creates a new matrix containing the gene of interest and
#' each binding partner that it interacts with.
#' @param MAE Input MAE which stores results from multiReg It is
#' recommended to use the MAE which was used in matrixFilter.
#' @param gene_interest Name of gene of interest. Full string required with no
#' spelling mistakes e.g. "mmu-miR-140-5p", "Ffg1", "hsa-miR-29a-5p". If gene
#' of interest is an mRNA, mRNAreg must be TRUE. Otherwise, if the gene of interest
#' is a miRNA, then mRNAreg must be FALSE.
#' @param mRNAreg TRUE or FALSE. Is the gene of interest a mRNA? Default is TRUE.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data
#' @return A matrix which contains the gene of interest and all binding partners.
#' Their values (Log2FC or ave exp) for each time point are also produced.
#' @export
#' @usage multiReg(MAE, gene_interest, mRNAreg, filt_df, miRNA_exp, mRNA_exp)
#' @importFrom stringr str_detect
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
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus', )
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
multiReg <- function(MAE, gene_interest, mRNAreg = TRUE,
                     filt_df, miRNA_exp, mRNA_exp){

  if (missing(MAE)) stop('MAE is missing. Add MAE to store output of multiReg.')

  if (missing(gene_interest)) stop('gene_interest is missing. Add string of gene of interest.')

  if (missing(filt_df)) stop('filt_df is missing. Add assay/ dataframe created by the matrifFilter function.')

  if (missing(miRNA_exp)) stop('miRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on miRNA expression data/ DE data.')

  if (missing(mRNA_exp)) stop('mRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on mRNA expression data/ DE data.')

  mRNA <- miR <- NULL

  x <- miRNA_exp

  x$ID <- NULL

  if(length(colnames(x)) < 5) {
    print('Warning: Fewer than five time points detected. This dataset is not suitable for regression based predictions! Results may be overestimated.')
  }

  Data <- rbind(miRNA_exp, mRNA_exp)

  Data$ID <- NULL

  Interactions <- filt_df

  Data1 <- Data[which(rownames(Data) %in% Interactions$mRNA),]

  Data2 <- Data[which(rownames(Data) %in% Interactions$miR),]

  if (mRNAreg == TRUE) {

    Data <- rbind(Data1, Data2)

    Data <- as.matrix(Data)

    Ints <- Interactions %>% filter(str_detect(mRNA, gene_interest))

    Ints <- Ints[which(Ints$mRNA %in% gene_interest == TRUE),]

    Vec <- as.array(c(unique(Ints$mRNA), Ints$miR))

  }else if (mRNAreg == FALSE) {

    Data <- rbind(Data2, Data1)

    Data <- as.matrix(Data)

    Ints <- Interactions %>% filter(str_detect(miR, gene_interest))

    Ints <- Ints[which(Ints$miR %in% gene_interest == TRUE),]

    Vec <- as.array(c(unique(Ints$miR), Ints$mRNA))

  }

  DF <- Data[which(rownames(Data) %in% Vec == TRUE),]

  DF <- t(as.matrix(DF))

  Time <- as.integer(gsub(rownames(DF), pattern = "[^0-9.-]",

                          replacement = ""))
  Time <- as.numeric(Time)

  DF <- cbind(Time, DF)

  MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(gene_interest = DF))))

  MAE <- suppressWarnings(c(MAE, MAE2))

  return(MAE)

}
