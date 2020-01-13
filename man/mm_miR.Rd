\name{mm_miR}
\alias{mm_miR}
\docType{data}
\title{Mouse microRNA data}
\description{
Mouse Fibrosis microRNA data from differential expression.
Folic Acid was injected into mouse kidneys to induce fibrosis and
nephropathy, and measurements were taken prior to and 1, 2, 3, 7, 14
and 28 days after injection. Data was put through limma voom for analysis.
}
\usage{data("mm_miR")}
\format{
A data frame with 278 observations on the following 10 variables.
\describe{
\item{\code{D1.Log2FC}}{a numeric vector}
\item{\code{D1.adjPVal}}{a numeric vector}
\item{\code{D2.Log2FC}}{a numeric vector}
\item{\code{D2.adjPVal}}{a numeric vector}
\item{\code{D3.Log2FC}}{a numeric vector}
\item{\code{D3.adjPVal}}{a numeric vector}
\item{\code{D7.Log2FC}}{a numeric vector}
\item{\code{D7.adjPVal}}{a numeric vector}
\item{\code{D14.Log2FC}}{a numeric vector}
\item{\code{D14.adjPVal}}{a numeric vector}
}
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE613287}
}
\references{
M. Pavkovic, L. Pantano, C. V. Gerlach, S. Brutus, S. A. Boswell,
R. A. Everley, J. V. Shah, S. H. Sui, and V. S. Vaidya, “Multi omics a
nalysis of fibrotic kidneys in two mouse models,” Scientific data, vol. 6,
no. 1, p. 92, 2019. <https://www.nature.com/articles/s41597-019-0095-5>
}
\examples{
data(mm_miR)
## maybe str(mm_miR) ; plot(mm_miR) ...
}
\keyword{datasets}
