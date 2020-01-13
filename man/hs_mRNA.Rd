\name{hs_mRNA}
\alias{hs_mRNA}
\docType{data}
\title{Human mRNA data set}
\description{
Differential expression from Human breast cancer cells (MCF-7).
Normoxic conditions contrasted against 16H, 32H and 48H under hypoxia
mRNA data. Only 2000 genes have been kept in this file for speed and size.
}
\usage{data("hs_mRNA")}
\format{
A data frame with 2000 observations on the following 6 variables.
\describe{
\item{\code{H16.logFC}}{a numeric vector}
\item{\code{H16.adjPVal}}{a numeric vector}
\item{\code{H32.logFC}}{a numeric vector}
\item{\code{H32.adjPVal}}{a numeric vector}
\item{\code{H48.logFC}}{a numeric vector}
\item{\code{H48.adjPVal}}{a numeric vector}
}
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47534}
}
\references{
C. Camps, H. K. Saini, D. R. Mole, H. Choudhry, M. Reczko, J. A.
Guerra-Assunção, Y.-M.Tian, F. M. Buffa, A. L. Harris, A. G. Hatzigeorgiou,
et al., “Integrated analysis of microrna and mrna expression and
association with hif binding reveals the complexity of microrna expression
regulation under hypoxia,” Molecular cancer, vol. 13, no. 1, p. 28, 2014..
<https://molecular-cancer.biomedcentral.com/articles/10.1186/1476-4598-13-28>
}
\examples{
data(hs_mRNA)
## maybe str(hs_mRNA) ; plot(hs_mRNA) ...
}
\keyword{datasets}
