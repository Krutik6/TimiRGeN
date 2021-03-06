% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createClusters2.R
\name{createClusters2}
\alias{createClusters2}
\title{createClusters2}
\usage{
createClusters2(MAE, genetic_data, noClusters)
}
\arguments{
\item{MAE}{MultiAssayExperiment which will store the results from
createClusters2.}

\item{genetic_data}{A dataframe with miR and mRNA
information together. This is the output from the combineGenes function
and will be stored as an assay within the MAE used in the combineGenes
function.}

\item{noClusters}{How many clusters should be generated? Default is 5.}
}
\value{
3 new objects in the input MAE.
Clusters(metadata): A list to be used as the input in checkClusters
and quickFuzz.
MfuzzData(ExperimentList): An ExpressionSet object to be used as input for
quickFuzz.
ClusterData(assay): An assay to be used as input for returnCluster.
}
\description{
Creates clusters from a dataframe of mRNAs and miRs. This
function should primarily be used when analysing data that has not gone
through pair-wise DE. This function will create clusters from longitudinal
temporal patterns. createClusters2 will create 3 data files.
1) Clusters will contain cluster logistics information and will be stored as
metadata, 2) MfuzzData will contain fuzzy clustering information and will be
stored as an experiment, 3) ClusterData will contain cluster-pathway fit
information and will be stored as an assay.
}
\examples{
data(long_data)
miRNA <- long_data[c(1:105),]
mRNA <- long_data[-c(1:105),]

MAE <- startObject(miR = miRNA, mRNA = mRNA)

MAE <- combineGenes(MAE, miR_data = assay(MAE, 1),
                    mRNA_data = assay(MAE, 2))

MAE <- createClusters2(MAE = MAE, genetic_data = assay(MAE, 3))
}
