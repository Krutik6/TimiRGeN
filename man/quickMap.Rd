% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quickMap.R
\name{quickMap}
\alias{quickMap}
\title{quickMap}
\usage{
quickMap(filt_df, numpairs)
}
\arguments{
\item{filt_df}{Dataframe from the matrixFilter function.}

\item{numpairs}{Number of pairs to plot. Must be an integer more than 1.}
}
\value{
Heatmap of miRNA-mRNA pairs.
}
\description{
Generates a heatmap of all miRNA:mRNA binding pairs that
have been filtered. Pairs are ranked by decreasing correlation.
}
\examples{
Filt_df <- data.frame(row.names = c("mmu-miR-145a-3p:Adamts15",
                                   "mmu-miR-146a-5p:Acy1"),
                     corr = c(-0.9191653, -0.7826041),
                     miR = c("mmu-miR-145a-3p", "mmu-miR-146a-5p"),
                     mRNA = c("Adamts15", "Acy1"),
                     miR_Entrez = c(387163, NA),
                     mRNA_Entrez = c(235130, 109652),
                     TargetScan = c(1, 0),
                     miRDB = c(0, 0),
                     Predicted_Interactions = c(1, 0),
                     miRTarBase = c(0, 1),
                     Pred_Fun = c(1, 1))

MAE <- MultiAssayExperiment()

MAE <- matrixFilter(MAE, miningMatrix = Filt_df, negativeOnly = FALSE,
                   threshold = 1, predictedOnly = FALSE)

quickMap(filt_df = MAE[[1]], numpairs = 2)
}
