\name{mm_mRNA}
\alias{mm_mRNA}
\docType{data}
\title{Mouse microRNA data}
\description{
Mouse Fibrosis mRNA data from differential expression.
Folic Acid was injected into mouse kidneys to induce fibrosis and
nephropathy, and measurements were taken prior to and 1, 2, 3, 7 and 14
days after injection. The 0 time point was contrasted over other time points. Data was put through limma for analysis. mm_mRNA only has DE results from 2000 mRNA genes for speed and size optimisation.

Column names are in the following format: timepoint(days).DEresult (Log2FC or adjPVal).
}
\usage{data("mm_mRNA")}
\format{
A data frame of mRNAs put through differential expression with 2000 observations on the following 10 variables.
\describe{
\item{\code{D1.Log2FC}}{a numeric vector containing log2FC values of
D1_fibrosis/ D0_fibrosis.}
\item{\code{D1.adjPVal}}{a numeric vector containing adjusted P values of
D1_fibrosis/ D0_fibrosis.}
\item{\code{D2.Log2FC}}{a numeric vector containing log2FC values of
D2_fibrosis/ D0_fibrosis.}
\item{\code{D2.adjPVal}}{a numeric vector containing adjusted P values of
D2_fibrosis/ D0_fibrosis.}
\item{\code{D3.Log2FC}}{a numeric vector containing log2FC values of
D3_fibrosis/ D0_fibrosis.}
\item{\code{D3.adjPVal}}{a numeric vector containing adjusted P values of
D3_fibrosis/ D0_fibrosis.}
\item{\code{D7.Log2FC}}{a numeric vector containing log2FC values of
D7_fibrosis/ D0_fibrosis.}
\item{\code{D7.adjPVal}}{a numeric vector containing adjusted P values of
D7_fibrosis/ D0_fibrosis.}
\item{\code{D14.Log2FC}}{a numeric vector containing log2FC values of
D14_fibrosis/ D0_fibrosis.}
\item{\code{D14.adjPVal}}{a numeric vector containing adjusted P values of
D14_fibrosis/ D0_fibrosis.}
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
data(mm_mRNA)
## maybe str(mm_mRNA) ; plot(mm_mRNA) ...
}
\keyword{datasets}
