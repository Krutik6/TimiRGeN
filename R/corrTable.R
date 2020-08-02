#' @title corrTable
#' @description An internal function for the mirMranInt function. This will
#' work out the correlations between the time series data from
#' diffExpressRes functions.
#' @param names Dataframe with combined list of mRNA (single wikipathway of
#'  interest) and miR names.
#' @param miR_exprs miR averaged expression data.
#' @param mRNA_wiki Single wikipathway of interest's mRNA averaged expression
#' data.
#' @param maxInt Maximum number of data samples in your datasets e.g. how many
#' timepoints.
#' @return A matrix of each mRNA and miR interaction and a correlation.
#' @importFrom stats cor
#' @noRd
corrTable <- function(names, miR_exprs, mRNA_wiki, maxInt){

    if (missing(names)) stop('Input dataframe of all combinations of miRs and
                              mRNAs.')

    if (missing(miR_exprs)) stop('Please use the diffExpressRes function
                                 on miR data. Will
                                 be stored in a MAE object as an assay.')

    if (missing(mRNA_wiki)) stop('Please use the wikiMrna function. Will
                                 be stored in a MAE object as an assay.')

    if (missing(maxInt)) stop('Input interger. Number of samples in your data.'
                              )

    # Each miR and mRNA could interact so extract their names
    n1 <- as.character(names[1])

    n2 <- as.character(names[2])

    # Create average correlation
    r <- stats::cor(as.numeric(miR_exprs[n1, seq_len(maxInt)]),
                    as.numeric(mRNA_wiki[n2, seq_len(maxInt)]))

    # Create correlation matrix
    data.frame(row.names = paste(n1,":",n2, sep = ""),
               avecor =r,
               miR = names[1],
               mRNA = names[2],
               miR_Entrez = miR_exprs[n1, -c(seq_len(maxInt))],
               mRNA_Entrez = mRNA_wiki[n2, -c(seq_len(maxInt))])
}

