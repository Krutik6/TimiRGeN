% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffExpressRes.R
\name{diffExpressRes}
\alias{diffExpressRes}
\title{diffExpressRes}
\usage{
diffExpressRes(MAE, df, dataType = '', genes_ID, idColumn = '',
name = '')
}
\arguments{
\item{MAE}{MultiAssayExperiment to store the output of diffExpressRes
within it. This function is to be used after pathways of interest have been
identified by enrichWiki or returnCluster. It is recommended to
store all diffExpressRes results in the MAE used in enrichWiki and/ or
returnCluster.}

\item{df}{mRNA or miR dataframe (rownames as genes and DE results as
columns). These will be found as assays in the MAE object used within the
startObject function.}

\item{dataType}{Column name to take an average from e.g. "Log2FC", "AveExp".
This string should be found consistently in the column names of your input
data. It is recommended to use a DE result value which represents abundance,
rather than confidence.}

\item{genes_ID}{Dataframe that was created from a getIds function e.g.
mRNA_ensembl or miR_entrez. Use the same ID type for miR and mRNA data.
These dataframes will be found as assays within the MAE which stores results
from the getIds functions.}

\item{idColumn}{Name of column to use as the merge point. If Column names in
getIds results have not been changed, it should be "GENENAME". Default has
been left as "GENENAME".}

\item{name}{New name of the assay. Should be a unique string. Remember
each assay in a MAE must have a unique name.}
}
\value{
Dataframe with only a single result type from DE (e.g. Log2FC) and
an ID type e.g. entrezIDs. Output will be stored as an assay in the input
MAE.
}
\description{
diffExpressRes will produce a dataframe which contains data
for only one result type, along with an ID of choice. It is recommended to
use this function on a DE results which represents abundance such as log2fc
or average expression, as this data will be averaged and correlated later in
the analysis. This is to be used for miR and mRNA data individually.
}
\examples{
library(org.Mm.eg.db)

miR <- mm_miR[1:50,]

mRNA <- mm_mRNA[1:100,]

MAE <- startObject(miR = miR, mRNA = mRNA)

MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')

MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus', orgDB = org.Mm.eg.db)

MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
                     genes_ID = assay(MAE, 7),
                     idColumn = 'GENENAME',
                     name = "mRNA_log2fc")
}
