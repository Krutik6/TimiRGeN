\name{Hs_gmt}
\alias{Hs_gmt}
\docType{data}
\title{
Hs_gmt
}
\description{
Hs_gmt
}
\usage{data("Hs_gmt")}
\format{
The format is:
Formal class 'MultiAssayExperiment' [package "MultiAssayExperiment"] with
5 slots
..@ ExperimentList:Formal class 'ExperimentList' [package
"MultiAssayExperiment"] with 4 slots
.. .. ..@ listData       :List of 3
.. .. .. ..$ path_gene:'data.frame':	24712 obs. of  2 variables:
.. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
.. .. .. .. ..$ gene: chr [1:24712] "574413" "2167" "4690" "5781" ...
.. .. .. ..$ path_name:'data.frame':	24712 obs. of  2 variables:
.. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
.. .. .. .. ..$ name: chr [1:24712] "FABP4 in ovarian cancer" "FABP4 in 
ovarian cancer" "B Cell Receptor Signaling Pathway" "B Cell Receptor Signaling 
Pathway" ...
.. .. .. ..$ path_data:'data.frame':	24712 obs. of  3 variables:
.. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
.. .. .. .. ..$ gene: chr [1:24712] "574413" "2167" "4690" "5781" ...
.. .. .. .. ..$ name: chr [1:24712] "FABP4 in ovarian cancer" "FABP4 in
ovarian cancer" "B Cell Receptor Signaling Pathway" "B Cell Receptor Signaling
Pathway" ...
.. .. ..@ elementType    : chr "ANY"
.. .. ..@ elementMetadata: NULL
.. .. ..@ metadata       : list()
..@ colData       :Formal class 'DataFrame' [package "S4Vectors"] with 6
slots
.. .. ..@ rownames       : chr(0) 
.. .. ..@ nrows          : int 0
.. .. ..@ listData       : Named list()
.. .. ..@ elementType    : chr "ANY"
.. .. ..@ elementMetadata: NULL
.. .. ..@ metadata       : list()
..@ sampleMap     :Formal class 'DataFrame' [package "S4Vectors"] with 6 
slots
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
Human gmt data for use in the vignette.
}
\source{
"http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Homo_sapiens.gmt"
}
\references{
"http://data.wikipathways.org/20200110/gmt/",
"wikipathways-20200110-gmt-Homo_sapiens.gmt"
}
\examples{
data(Hs_gmt)
## maybe str(Hs_gmt) ; plot(Hs_gmt) ...
}
\keyword{datasets}
