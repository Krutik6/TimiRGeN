\name{long_data}
\alias{long_data}
\docType{data}
\title{
Data from a longer time series (10 time points)
}
\description{
Averaged raw count data of miRNAs (105) and mRNAs (1041). Time based differential expression was performed with DESeq2 (using LTR method). Differentially expressed genes were filtered from average raw counts. miRNA and mRNA data was combined into one dataframe. This is only example data used in section 5 of the
vignette.
}
\usage{data("long_data")}
\format{
  A data frame with 1146 observations on the following 10 variables.
  \describe{
    \item{\code{H0.Counts}}{a numeric vector of Raw counts at 0 hours.}
    \item{\code{H1.Counts}}{a numeric vector of Raw counts at 1 hours.}
    \item{\code{H2.Counts}}{a numeric vector of Raw counts at 2 hours.}
    \item{\code{H3.Counts}}{a numeric vector of Raw counts at 3 hours.}
    \item{\code{H4.Counts}}{a numeric vector of Raw counts at 4 hours.}
    \item{\code{H5.Counts}}{a numeric vector of Raw counts at 5 hours.}
    \item{\code{H6.Counts}}{a numeric vector of Raw counts at 6 hours.}
    \item{\code{H8.Counts}}{a numeric vector of Raw counts at 8 hours.}
    \item{\code{H12.Counts}}{a numeric vector of Raw counts at 12 hours.}
    \item{\code{H24.Counts}}{a numeric vector of Raw counts at 24 hours.}
  }
}
\references{
Baran-Gale J, Purvis JE, Sethupathy P. An integrative transcriptomics approach identifies miR-503 as a candidate master regulator of the estrogen response in MCF-7 breast cancer cells. RNA 2016 Oct;22(10):1592-603. PMID: 27539783
}
\examples{
data(long_data)
miRNA <- long_data[c(1:105),]
mRNA <- long_data[-c(1:105),]
}
\keyword{datasets}
