\name{w_list_mouse}
\alias{w_list_mouse}
\docType{data}
\title{
Wikipathways lists for mouse produced by the wikiList function.
}
\description{
List of entrezIDs associated with each wikipathway made for mouse. This dataset
will allow for functions to run faster during CMD check and vignette building.
This is used instead of the wikiList function in several instances in TimiRGeN
because the wikiList takes a lot of time to download data. For example in the
vignette and in several function examples. To reproduce the output enter the
following code.
>MAE <- MultiAssayExperiment
>MAE2 <- wikiList(MAE, stringSpecies = 'Mus musculus', stringSymbol = 'L')
}
\usage{data("w_list_mouse")}
\format{
List of every mouse wikipathway. Each wikipathway will have associated genes
stored as entrezIDs.
}
\details{
A Large list of each wikipathway and associated enterzgene IDs.
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
data(w_list_mouse) -> w_list_mouse
}
\keyword{datasets}
