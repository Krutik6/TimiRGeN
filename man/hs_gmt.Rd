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
 A MultiAssayExperiment object of 5 listed
 experiments with user-defined names and respective classes. 
 Containing an ExperimentList class object of length 5: 
 [1] path_gene: data.frame with 24712 rows and 2 columns 
 [2] path_name: data.frame with 24712 rows and 2 columns 
 [3] path_data: data.frame with 24712 rows and 3 columns 
 [4] path_gene_ensembl: data.frame with 27466 rows and 2 columns 
 [5] path_data_ensembl: data.frame with 27466 rows and 3 columns 
Features: 
 experiments() - obtain the ExperimentList instance 
 colData() - the primary/phenotype DataFrame 
 sampleMap() - the sample availability DataFrame 
 `$`, `[`, `[[` - extract colData columns, subset, or experiment 
 *Format() - convert into a long or wide DataFrame 
 assays() - convert ExperimentList to a SimpleList of matrices
 metadata()
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
