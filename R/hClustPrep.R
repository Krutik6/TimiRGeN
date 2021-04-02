#' @title hClustPrep
#' @description Internal function for further preparation for
#' hierarchical clustering.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param distmeth Dist method for hierarchical clustering. Default is
#' "maximum". Type >help(dist) in console for more options.
#' @param hclustmeth Hclust method for hierarchical clustering. Default is
#' "ward.D". Type >help(hclust) in console for more options.
#' @noRd
#' @importFrom stats dist hclust
#' @importFrom tidyr drop_na spread
hClustPrep <- function(filt_df, miRNA_exp, mRNA_exp, distmeth="maximum",
                       hclustmeth = "ward.D"){

  if (missing(filt_df)) stop('filt_df is missing. Add assay/ dataframe created by the matrixFilter function.')

  if (missing(miRNA_exp)) stop('miRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on miRNA expression data/ DE data.')

  if (missing(mRNA_exp)) stop('mRNA_exp is missing Add assay/ dataframe created by the diffExpressRes function used on mRNA expression data/ DE data.')

  spread <- Gene <- Expression <- NULL

  Prep <- clustPrep(filt_df, miRNA_exp, mRNA_exp)

  NaD <- drop_na(Prep)

  spread_Genes <- NaD %>% spread(Gene, Expression)

  newSpread <- t(spread_Genes[-1])

  dists <- stats::dist(newSpread, method=distmeth)

  fit <- stats::hclust(dists, method=hclustmeth)

  return(fit)

}
