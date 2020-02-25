\name{Mm_gmt}
\alias{Mm_gmt}
\docType{data}
\title{
Mm_gmt
}
\description{
mouse gmt data
}
\usage{data("Mm_gmt")}
\format{
The format is:
Formal class 'MultiAssayExperiment' [package "MultiAssayExperiment"] with 
5 slots
..@ ExperimentList:Formal class 'ExperimentList' [package
"MultiAssayExperiment"] with 4 slots
.. .. ..@ listData       :List of 3
.. .. .. ..$ path_gene:'data.frame':	9786 obs. of  2 variables:
.. .. .. .. ..$ wpid: chr [1:9786] "WP571" "WP571" "WP571" "WP571" ...
.. .. .. .. ..$ gene: chr [1:9786] "16175" "12370" "26419" "19645" ...
.. .. .. ..$ path_name:'data.frame':	9786 obs. of  2 variables:
.. .. .. .. ..$ wpid: chr [1:9786] "WP571" "WP571" "WP571" "WP571" ...
.. .. .. .. ..$ name: chr [1:9786] "FAS pathway and Stress induction of HSP
regulation" "FAS pathway and Stress induction of HSP regulation" "FAS pathway 
and Stress induction of HSP regulation" "FAS pathway and Stress induction of 
HSP regulation" ...
.. .. .. ..$ path_data:'data.frame':	9786 obs. of  3 variables:
.. .. .. .. ..$ wpid: chr [1:9786] "WP571" "WP571" "WP571" "WP571" ...
.. .. .. .. ..$ gene: chr [1:9786] "16175" "12370" "26419" "19645" ...
.. .. .. .. ..$ name: chr [1:9786] "FAS pathway and Stress induction of HSP
regulation" "FAS pathway and Stress induction of HSP regulation" "FAS pathway
and Stress induction of HSP regulation" "FAS pathway and Stress induction of
HSP regulation" ...
.. .. ..@ elementType    : chr "ANY"
.. .. ..@ elementMetadata: NULL
.. .. ..@ metadata       : list()
..@ colData       :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
.. .. ..@ rownames       : chr(0) 
.. .. ..@ nrows          : int 0
.. .. ..@ listData       : Named list()
.. .. ..@ elementType    : chr "ANY"
.. .. ..@ elementMetadata: NULL
.. .. ..@ metadata       : list()
..@ sampleMap     :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
.. .. ..@ rownames       : NULL
.. .. ..@ nrows          : int 0
.. .. ..@ listData       :List of 3
.. .. .. ..$ assay  : Factor w/ 0 levels: 
.. .. .. ..$ primary: chr(0) 
.. .. .. ..$ colname: chr(0) 
.. .. ..@ elementType    : chr "ANY"
.. .. ..@ elementMetadata: NULL
.. .. ..@ metadata       : list()
..@ drops         : list()
..@ metadata      : list()
}
\details{
Mouse gmt data availabel for use in the vignette.
}
\source{
"http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Mus_musculus.gmt"
}
\references{
"http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Mus_musculus.gmt"
}
\examples{
data(Mm_gmt)
## maybe str(Mm_gmt) ; plot(Mm_gmt) ...
}
\keyword{datasets}
