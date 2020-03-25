\name{mmu_miRTarBase}
\alias{mmu_miRTarBase}
\docType{data}
\title{
mmu_miRTarBase
}
\description{
mmu_miRTarBase. Used instead of downloading data while their servers are down.
}
\usage{data("mmu_miRTarBase")}
\format{
A data frame with 1629 observations on the following 3 variables.
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
data(mmu_miRTarBase)
## maybe str(mmu_miRTarBase) ; plot(mmu_miRTarBase) ...
}
\keyword{datasets}
