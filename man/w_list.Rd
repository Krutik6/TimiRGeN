\name{w_list}
\alias{w_list}
\docType{data}
\title{
Wikipathways lists for mouse
}
\description{
List of entrezIDs associated with each wikipathway made for mouse. This dataset
will allow for functions to run faster during CMD check and vignette building.
}
\usage{data("w_list")}
\format{
list
}
\details{
A Large list of each wikipathway and associated genes.
}
\source{
https://www.wikipathways.org/index.php/Download_Pathways
}
\references{
Denise N Slenter et al. “WikiPathways: a multifaceted pathway database b
ridging metabolomics to other omics research”.
In: Nucleic acids research 46.D1 (2017), pp. D661–D667.
}
\examples{
data(w_list) -> w_list
}
\keyword{datasets}
