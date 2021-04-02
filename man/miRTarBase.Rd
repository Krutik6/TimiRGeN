\name{miRTarBase}
\alias{miRTarBase}
\docType{data}
\title{
miRTarBase dataset
}
\description{
miRTarBase dataset - downloaded from http://mirtarbase.cuhk.edu.cn/php/download.php on 7/03/21. miRNA-mRNA interactions labelled as "weak" have been removed.
}
\usage{data("miRTarBase")}
\format{
  A data frame with 13315 observations on the following 9 variables.
  \describe{
    \item{\code{miRTarBase.ID}}{Column of characters describing the ID of each
    miRNA-mRNA interaction.}
    \item{\code{miRNA}}{Column of characters describing the ID of each
    miRNA from an interaction.}
    \item{\code{Species..miRNA.}}{Column of characters describing the name of
    the species from each miRNA of an interaction.}
    \item{\code{Target.Gene}}{Column of characters describing the ID of each
    mRNA from an interaction.}
    \item{\code{Target.Gene..Entrez.ID.}}{Column of characters describing the
    entrezgene ID of each miRNA from an interaction.}
    \item{\code{Species..Target.Gene.}}{Column of characters describing the name
    of the species from each mRNA of an interaction.}
    \item{\code{Experiments}}{Column of characters describing the experiment        used  to capture each miNRA-mRNA interaction.}
    \item{\code{Support.Type}}{Column of characters describing the strength of
    each experiment.}
    \item{\code{References..PMID.}}{Column of characters describing the         publication ID from each miRNA-mRNA interaction.}
  }
}
\source{
http://mirtarbase.cuhk.edu.cn/php/download.php
}
\references{
Huang HY, Lin YC, Li J, Huang KY, Shrestha S, Hong HC, Tang Y, Chen YG, Jin CN, Yu Y, Xu JT, Li YM, Cai XX, Zhou ZY, Chen XH, Pei YY, Hu L, Su JJ, Cui SD, Wang F, Xie YY, Ding SY, Luo MF, Chou CH, Chang NW, Chen KW, Cheng YH, Wan XH, Hsu WL, Lee TY, Wei FX, Huang HD*   "miRTarBase 2020: updates to the experimentally validated microRNA-target interaction database"  Nucleic Acids Research 2020 Jan 8;48(D1):D148-D154.
}
\examples{
data(miRTarBase)
}
\keyword{datasets}
