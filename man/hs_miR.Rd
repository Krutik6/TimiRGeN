\name{hs_miR}
\alias{hs_miR}
\docType{data}
\title{Human microRNA data set}
\description{
Differential expression from Human breast cancer cells (MCF-7).
Normoxic conditions contrasted against 16H, 32H and 48H under hypoxia. Data was put though limma for DE. hs_miR consists DE results from 189 microRNAs.

Column names are in the following format: timepoint(hours).DEresult (logFC or adjPVal).

The miR names are not using in TimiRGeN friendly nomenclature. Gene names must be changed before use in TimiRGeN analysis. This must be changed before using TimiRGeN. miR gene name changing is explained in section 3.1 of the vignette.
}
\usage{data("hs_miR")}
\format{
A data frame of miRNAs put through differential expression with 189 observations on the following 6 variables.
\describe{
\item{\code{H16.logFC}}{a numeric vector containing log2FC values of
Hypoxia_16_hours/ Normoxic conditions.}
\item{\code{H16.adjPVal}}{a numeric vector containing adjusted P values values
of Hypoxia_16_hours/ Normoxic conditions.}
\item{\code{H32.logFC}}{a numeric vector containing log2FC values of
Hypoxia_32_hours/ Normoxic conditions.}
\item{\code{H32.adjPVal}}{a numeric vector containing adjusted P values values
of Hypoxia_32_hours/ Normoxic conditions}
\item{\code{H48.logFC}}{a numeric vector containing log2FC values of
Hypoxia_48_hours/ Normoxic conditions.}
\item{\code{H48.adjPVal}}{a numeric vector containing adjusted P values values
of Hypoxia_48_hours/ Normoxic conditions}
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
data(hs_miR)
## maybe str(hs_miR) ; plot(hs_miR) ...
}
\keyword{datasets}
