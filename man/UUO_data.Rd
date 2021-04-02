\name{UUO_data}
\alias{UUO_data}
\docType{data}
\title{
UUO_data
}
\description{
An alternative mouse kidney fibrosis dataset which is used in section 6 of the
vignette to demonstrate meta-analysis with TimiRGeN. UUO_data contains 240
miRNAs and 1874 mRNAs, all of which can also be found in mm_mRNA and mm_miR
data files. This subset was loaded into the package for speed. Log2FC and
adjusted P values for each gene is found. Differential expression was conducted
on a pair-wise manor (D3/D0, D7/D0 and D14/D0). Unlike the folic acid data
(mm_miR and mm_mRNA), the UUO data did not include time points D1 and D2.

}
\usage{data("UUO_data")}
\format{
  A data frame with 2114 observations on the following 6 variables.
  \describe{
    \item{\code{D3.Log2FC}}{a numeric vector of log2fc values from D3/D0.
    Contains miRNA and mRNA data.}
    \item{\code{D3.ajdPVal}}{a numeric vector of adj p values values from D3/D0.
    Contains miRNA and mRNA data.}
    \item{\code{D7.Log2FC}}{a numeric vector of log2fc values from D7/D0.
    Contains miRNA and mRNA data.}
    \item{\code{D7.ajdPVal}}{a numeric vector of adj p values values from D7/D0.
    Contains miRNA and mRNA data.}
    \item{\code{D14.Log2FC}}{a numeric vector of log2fc values from D14/D0.
    Contains miRNA and mRNA data.}
    \item{\code{D14.ajdPVal}}{a numeric vector of adj p values values from
    D14/D0. Contains miRNA and mRNA data.}
  }
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118341
}
\references{
Pavkovic M, Pantano L, Gerlach CV, Brutus S et al. Multi omics analysis of fibrotic kidneys in two mouse models. Sci Data 2019 Jun 14;6(1):92.
}
\examples{
data(UUO_data)
UUO_miRNA <- UUO_data[c(1:240),]
UUO_mRNA <- UUO_data[-c(1:240),]
}
\keyword{datasets}
