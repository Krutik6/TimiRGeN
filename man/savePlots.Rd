% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/savePlots.R
\name{savePlots}
\alias{savePlots}
\title{savePlots}
\usage{
savePlots(largeList, maxInt, fileType = '', width, height)
}
\arguments{
\item{largeList}{A large list containing GSEA results. This should be
stored as metadata within the MAE used in the enrichWiki function.}

\item{maxInt}{Integer, number of samples in data set.}

\item{fileType}{Type of file for images to be exported as: "png", "tiff",
"svg" or "jpeg".}

\item{width}{= Width of plots in inches. Default is 22 inches.}

\item{height}{= Heightof plots in inches. Default is 10 inches.}
}
\value{
Saves plots in working directory. Each sample (e.g. time point) will
have a separate plot.
}
\description{
Saves all plots from enrichWiki into the current working
directory.
}
