#' @title corrTable
#' @description An internal function for the mirMrnaInt function. This will
#' work out the pearson correlations between the miR and mRNA DE results data
#' over the time series.
#' @param names Dataframe with combined list of mRNA (single wikipathway of
#'  interest) and miR names. Produced by mirMrnaInt.
#' @param miR_exprs miR dataframe with a single DE result type (e.g. log2fc) and
#' an ID type.
#' @param mRNA_wiki mRNA dataframe containing mRNAs found in the pathway
#' of interest. Will also contain a single DE result type (e.g. log2fc) and
#' an ID type.
#' @param maxInt Maximum number of data samples in your datasets e.g. how many
#' time points.
#' @return A matrix of each mRNA and miR interaction and averaged correlations.
#' @importFrom stats cor
#' @noRd
corrTable <- function(names, miR_exprs, mRNA_wiki, maxInt){

    # Each miR and mRNA could interact so extract their names
    n1 <- as.character(names[1])

    n2 <- as.character(names[2])

    # Create correlation
    r <- stats::cor(as.numeric(miR_exprs[n1, seq_len(maxInt)]),
                    as.numeric(mRNA_wiki[n2, seq_len(maxInt)]))

    # Create correlation matrix
    data.frame(row.names = paste(n1,":",n2, sep = ""),
               corr =r,
               miR = names[1],
               mRNA = names[2],
               miR_Entrez = miR_exprs[n1, -c(seq_len(maxInt))],
               mRNA_Entrez = mRNA_wiki[n2, -c(seq_len(maxInt))])
}

