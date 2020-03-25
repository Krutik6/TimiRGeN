\name{Int_data}
\alias{Int_data}
\docType{data}
\title{
Int_data
}
\description{
Data frame of the correlations found between the DE lof2fc scores from the
miR and mRNA data from Mouse Kidney Fibrosis data set. Includes if these
interactions were also found in TargetScans, miRDB and miRTarBase.
}
\usage{data("Int_data")}
\format{
A data frame with 6394 observations on the following 10 variables.
\describe{
\item{\code{avecor}}{a numeric vector}
\item{\code{miR}}{a factor vector}
\item{\code{mRNA}}{a factor vector}
\item{\code{miR_Entrez}}{a factor vector}
\item{\code{mRNA_Entrez}}{a factor vector}
\item{\code{TargetScan}}{a numeric vector}
\item{\code{miRDB}}{a numeric vector}
\item{\code{Predicted_Interactions}}{a numeric vector}
\item{\code{miRTarBase}}{a numeric vector}
\item{\code{Pred_Fun}}{a numeric vector}
}
}
\details{
Dataframe of matrix_mining after DataMiningMatrix function. To use for 
vignette building.
}
\source{
TimiRGeN package.
}
\references{
NA
}
\examples{
data(Int_data)
## maybe str(Int_data) ; plot(Int_data) ...
}
\keyword{datasets}
