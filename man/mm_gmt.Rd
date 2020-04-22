\name{mm_gmt}
\alias{mm_gmt}
\docType{data}
\title{
mm_gmt
}
\description{
mouse gmt data
}
\usage{data("mm_gmt")}
\format{
A MultiAssayExperiment object of 3 listed
 experiments with user-defined names and respective classes. 
 Containing an ExperimentList class object of length 3: 
 [1] path_gene: data.frame with 9786 rows and 2 columns 
 [2] path_name: data.frame with 9786 rows and 2 columns 
 [3] path_data: data.frame with 9786 rows and 3 columns 
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
data(mm_gmt)
## maybe str(mm_gmt) ; plot(mm_gmt) ...
}
\keyword{datasets}
