\name{e_list_mouse}
\alias{e_list_mouse}
\docType{data}
\title{
e_list_mouse
}
\description{
List of entrezIDs of the significantly DE genes from the mouse fibrosis dataset. To make examples run faster this data is used, instead of re-running previous functions
throughout the examples of TimiRGeN. e_list_mouse is the output of the eNames function
when the combined analysis is performed on mm_miR and mm_RNA and when entrezIDs
are used as the gene IDs.
}
\usage{data("e_list_mouse")}
\format{
The format is:
List of 5:
D1
D2
D3
D7
D14
Each list contains significantly differentially expressed entrezgene IDs,
specific for each time point.
}
\details{
List of entrezgeneIDs per time point (5) from combined analysis of mm_miR and
mm_mRNA. Used to speed up examples for building and checks.
}
\source{
From using eNames during combined analysis of mm_miR and mm_mRNA..
}
\references{
NA
}
\examples{
data(e_list_mouse)
}
\keyword{datasets}
