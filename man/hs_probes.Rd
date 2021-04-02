\name{hs_probes}
\alias{hs_probes}
\docType{data}
\title{hs_probes}
\description{
Gene names and ensemble genes IDs stored as characters. This is used as a demonstration in using a "user made" universe for overenrichment analysis. This is only used in section 3 of the vignette.}
\usage{data("hs_probes")}
\format{
  A data frame with 16801 observations on the following 2 variables.
  \describe{
    \item{\code{external_gene_name}}{a character vector of gene names of probes.}
    \item{\code{ensembl_gene_id}}{a character vector of ensemble gene IDs of
    probes.}
  }
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8227
}
\references{
Camps C, Saini HK, Mole DR, Choudhry H et al. Integrated analysis of microRNA and mRNA expression and association with HIF binding reveals the complexity of microRNA expression regulation under hypoxia. Mol Cancer 2014 Feb 11;13:28. PMID: 24517586
}
\examples{
data(hs_probes)
}
\keyword{datasets}
