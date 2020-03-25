\name{hs_gmt}
\alias{hs_gmt}
\docType{data}
\title{
hs_gmt
}
\description{
hs_gmt
}
\usage{data("hs_gmt")}
\format{
  The format is:
Formal class 'MultiAssayExperiment' [package "MultiAssayExperiment"] with 5 slots
  ..@ ExperimentList:Formal class 'ExperimentList' [package "MultiAssayExperiment"] with 4 slots
  .. .. ..@ listData       :List of 5
  .. .. .. ..$ path_gene        :'data.frame':	24712 obs. of  2 variables:
  .. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
  .. .. .. .. ..$ gene: chr [1:24712] "574413" "2167" "4690" "5781" ...
  .. .. .. ..$ path_name        :'data.frame':	24712 obs. of  2 variables:
  .. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
  .. .. .. .. ..$ name: chr [1:24712] "FABP4 in ovarian cancer" "FABP4 in ovarian cancer" "B Cell Receptor Signaling Pathway" "B Cell Receptor Signaling Pathway" ...
  .. .. .. ..$ path_data        :'data.frame':	24712 obs. of  3 variables:
  .. .. .. .. ..$ wpid: chr [1:24712] "WP4400" "WP4400" "WP23" "WP23" ...
  .. .. .. .. ..$ gene: chr [1:24712] "574413" "2167" "4690" "5781" ...
  .. .. .. .. ..$ name: chr [1:24712] "FABP4 in ovarian cancer" "FABP4 in ovarian cancer" "B Cell Receptor Signaling Pathway" "B Cell Receptor Signaling Pathway" ...
  .. .. .. ..$ path_gene_ensembl:'data.frame':	27466 obs. of  2 variables:
  .. .. .. .. ..$ wpid: chr [1:27466] "WP694" "WP3633" "WP702" "WP3591" ...
  .. .. .. .. ..$ gene: chr [1:27466] "ENSG00000156006" "ENSG00000156006" "ENSG00000156006" "ENSG00000196839" ...
  .. .. .. ..$ path_data_ensembl:'data.frame':	27466 obs. of  3 variables:
  .. .. .. .. ..$ wpid.x: chr [1:27466] "WP694" "WP3633" "WP702" "WP3591" ...
  .. .. .. .. ..$ gene  : chr [1:27466] "ENSG00000156006" "ENSG00000156006" "ENSG00000156006" "ENSG00000196839" ...
  .. .. .. .. ..$ name  : chr [1:27466] "Arylamine metabolism" "Caffeine and Theobromine metabolism" "Metapathway biotransformation Phase I and II" "Sleep regulation" ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ colData       :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : chr [1:4] "wpid" "gene" "name" "wpid.x"
  .. .. ..@ nrows          : int 4
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ sampleMap     :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 12
  .. .. ..@ listData       :List of 3
  .. .. .. ..$ assay  : Factor w/ 5 levels "path_gene","path_name",..: 1 1 2 2 3 3 3 4 4 5 ...
  .. .. .. ..$ primary: chr [1:12] "wpid" "gene" "wpid" "name" ...
  .. .. .. ..$ colname: chr [1:12] "wpid" "gene" "wpid" "name" ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ drops         : list()
  ..@ metadata      : list()
}
\details{
Human gmt data availabel for use in the vignette.
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
data(hs_gmt)
## maybe str(hs_gmt) ; plot(hs_gmt) ...
}
\keyword{datasets}
