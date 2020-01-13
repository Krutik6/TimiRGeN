#' @title CorrTable
#'
#' @param names Dataframe with combined list of mRNA(single wikipathway of
#'  interest) and miR names.
#' @param miR_exprs miR averaged expression data.
#' @param mRNA_wiki Single wikipathway of interest's mRNA averaged expression
#' data.
#' @param maxInt Maximum number of data samples in your datasets e.g. how many
#'  timepoints.
#' @return A matrix of each mRNA and miR interaction and a correlation.
#' @export
#' @importFrom stats ave complete.cases median na.omit p.adjust reorder
CorrTable <- function(names, miR_exprs, mRNA_wiki, maxInt){
if (missing(names)) stop('Input dataframe of all combinations of miRs and
mRNAs.');
if (missing(miR_exprs)) stop('Input dataframe of miR data from the Express
function.');
if (missing(mRNA_wiki)) stop('Input dataframe of mRNA data found in
wikipathway of choice.');
if (missing(maxInt)) stop('Input interger. Number of samples in your data.');
n1 <- as.character(names[1])
n2 <- as.character(names[2])
r <- stats::cor(as.numeric(miR_exprs[n1, seq_len(maxInt)]),
as.numeric(mRNA_wiki[n2, seq_len(maxInt)]))
data.frame(row.names = paste(n1,":",n2, sep = ""),avecor =r,
miR = names[1], mRNA = names[2],
miR_Entrez = miR_exprs[n1, -c(seq_len(maxInt))],
mRNA_Entrez = mRNA_wiki[n2, -c(seq_len(maxInt))])
}

