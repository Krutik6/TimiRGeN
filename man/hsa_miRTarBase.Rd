\name{hsa_miRTarBase}
\alias{hsa_miRTarBase}
\docType{data}
\title{
hsa_miRTarBase
}
\description{
hsa_miRTarBase. Used instead of downloading data while their servers are down.
}
\usage{data("hsa_miRTarBase")}
\format{
  A data frame with 10754 observations on the following 3 variables.
\describe{
\item{\code{miRTarBase_Interactions}}{a factor vector}
\item{\code{miRTarBase_microRNA}}{a factor vector}
\item{\code{miRTarBase_mRNA}}{a factor vector}
}
}
\details{
miRTarBase evidence that is not labelled "weak" for human miR-mRNA interactions.
}
\source{
http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/miRTarBase_MTI.xlsx
}
\references{
http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/miRTarBase_MTI.xlsx
}
\examples{
data(hsa_miRTarBase)
## maybe str(hsa_miRTarBase) ; plot(hsa_miRTarBase) ...
}
\keyword{datasets}
