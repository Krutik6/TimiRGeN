#' @title clustPrep
#' @description Internal function to prepare data for hierarchical clustering.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data.
#' @noRd
#' @importFrom reshape2 melt
clustPrep <- function(filt_df, miRNA_exp, mRNA_exp){

  if (missing(filt_df)) stop('MAE is missing. Add assay/ dataframe created by the matrixFilter function.')

  if (missing(miRNA_exp)) stop('miRNA_exp is missing Add assay/ dataframe created by the diffExpressRes function used on miRNA expression data/ DE data.')

  if (missing(mRNA_exp)) stop('mRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on mRNA expression data/ DE data.')

  Ranks <- filt_df[,c(1,2,3)][order(filt_df$corr, decreasing = FALSE),]

  miR_genes <- miRNA_exp[which(rownames(miRNA_exp) %in% Ranks$miR),]

  mRNA_genes <- mRNA_exp[which(rownames(mRNA_exp) %in% Ranks$mRNA),]

  genes <- rbind(miR_genes, mRNA_genes)

  genes$ID <- NULL

  genes <- scale(genes)

  colnames(genes) <- as.integer(gsub(colnames(genes), pattern = "[^0-9.-]",
                                     replacement = ""))

  hcprep <- reshape2::melt(as.matrix(genes), varnames = c("Gene", "Time"))

  names(hcprep)[[3]] <- "Expression"

  return(hcprep)

}
