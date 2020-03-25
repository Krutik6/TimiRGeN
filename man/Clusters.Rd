\name{Clusters}
\alias{Clusters}
\docType{data}
\title{
Clusters
}
\description{
Results from CreateClusters. This will speeden the process of vignette building
and avoid biocondoctor error.
}
\usage{data("Clusters")}
\format{
  The format is:
Formal class 'MultiAssayExperiment' [package "MultiAssayExperiment"] with 5 slots
  ..@ ExperimentList:Formal class 'ExperimentList' [package "MultiAssayExperiment"] with 4 slots
  .. .. ..@ listData       :List of 2
  .. .. .. ..$ ClusterData:'data.frame':	200 obs. of  13 variables:
  .. .. .. .. ..$ 1          : num [1:200] 0.08271 0.00246 0.07673 0.04898 0.05954 ...
  .. .. .. .. ..$ 2          : num [1:200] 0.15848 0.00373 0.07328 0.11898 0.15695 ...
  .. .. .. .. ..$ 3          : num [1:200] 0.047978 0.000748 0.058073 0.017995 0.033526 ...
  .. .. .. .. ..$ 4          : num [1:200] 0.0838 0.9718 0.1736 0.1634 0.1244 ...
  .. .. .. .. ..$ 5          : num [1:200] 0.055685 0.000631 0.032666 0.014194 0.026322 ...
  .. .. .. .. ..$ 6          : num [1:200] 0.0648 0.000737 0.033712 0.01966 0.024815 ...
  .. .. .. .. ..$ 7          : num [1:200] 0.08376 0.00267 0.0837 0.04363 0.08074 ...
  .. .. .. .. ..$ 8          : num [1:200] 0.100785 0.000618 0.032122 0.017053 0.037406 ...
  .. .. .. .. ..$ 9          : num [1:200] 0.049188 0.000918 0.050536 0.025298 0.024316 ...
  .. .. .. .. ..$ 10         : num [1:200] 0.06486 0.00454 0.16612 0.06639 0.05861 ...
  .. .. .. .. ..$ 11         : num [1:200] 0.10709 0.00299 0.1042 0.06772 0.2738 ...
  .. .. .. .. ..$ 12         : num [1:200] 0.10083 0.00814 0.11524 0.39665 0.0996 ...
  .. .. .. .. ..$ Description: chr [1:200] "Statin Pathway" "IL-9 Signaling Pathway" "Cholesterol Biosynthesis" "Selenium metabolism/Selenoproteins" ...
  .. .. .. ..$ MfuzzData  :Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  .. .. .. .. .. ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. .. .. .. .. .. ..@ name             : chr ""
  .. .. .. .. .. .. .. ..@ lab              : chr ""
  .. .. .. .. .. .. .. ..@ contact          : chr ""
  .. .. .. .. .. .. .. ..@ title            : chr ""
  .. .. .. .. .. .. .. ..@ abstract         : chr ""
  .. .. .. .. .. .. .. ..@ url              : chr ""
  .. .. .. .. .. .. .. ..@ pubMedIds        : chr ""
  .. .. .. .. .. .. .. ..@ samples          : list()
  .. .. .. .. .. .. .. ..@ hybridizations   : list()
  .. .. .. .. .. .. .. ..@ normControls     : list()
  .. .. .. .. .. .. .. ..@ preprocessing    : list()
  .. .. .. .. .. .. .. ..@ other            : list()
  .. .. .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. .. .. .. ..@ assayData        :<environment: 0x2e39c418> 
  .. .. .. .. .. ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. .. .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. .. .. .. .. .. ..$ labelDescription: chr(0) 
  .. .. .. .. .. .. .. ..@ data             :'data.frame':	5 obs. of  0 variables
  .. .. .. .. .. .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. .. .. .. ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. .. .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. .. .. .. .. .. ..$ labelDescription: chr(0) 
  .. .. .. .. .. .. .. ..@ data             :'data.frame':	200 obs. of  0 variables
  .. .. .. .. .. .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. .. .. .. ..@ annotation       : chr(0) 
  .. .. .. .. .. ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. .. .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. .. .. .. .. .. ..$ labelDescription: chr(0) 
  .. .. .. .. .. .. .. ..@ data             :'data.frame':	5 obs. of  0 variables
  .. .. .. .. .. .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. .. ..@ .Data:List of 4
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 3 6 0
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 2 44 0
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ colData       :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : chr [1:18] "1" "2" "3" "4" ...
  .. .. ..@ nrows          : int 18
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ sampleMap     :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 18
  .. .. ..@ listData       :List of 3
  .. .. .. ..$ assay  : Factor w/ 2 levels "ClusterData",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ primary: chr [1:18] "1" "2" "3" "4" ...
  .. .. .. ..$ colname: chr [1:18] "1" "2" "3" "4" ...
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ drops         : list()
  ..@ metadata      :List of 1
  .. ..$ Clusters:List of 7
  .. .. ..$ centers    : num [1:12, 1:5] 0.556 0.562 -0.696 0.442 -0.413 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:12] "1" "2" "3" "4" ...
  .. .. .. .. ..$ : chr [1:5] "D1" "D2" "D3" "D7" ...
  .. .. ..$ size       : int [1:12] 18 22 10 27 6 11 16 17 17 17 ...
  .. .. ..$ cluster    : Named int [1:200] 2 4 4 12 11 5 9 11 9 9 ...
  .. .. .. ..- attr(*, "names")= chr [1:200] "WP1" "WP10" "WP103" "WP108" ...
  .. .. ..$ membership : num [1:200, 1:12] 0.08271 0.00246 0.07673 0.04898 0.05954 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:200] "WP1" "WP10" "WP103" "WP108" ...
  .. .. .. .. ..$ : chr [1:12] "1" "2" "3" "4" ...
  .. .. ..$ iter       : num 100
  .. .. ..$ withinerror: num 0.0522
  .. .. ..$ call       : language cmeans(x = exprs(eset), centers = centers, method = "cmeans", m = m)
  .. .. ..- attr(*, "class")= chr "fclust"
}
\details{
Bioconductor vignette build produces an unkown error, to get around this 
the following data will be added instead. This will also hasten the process.
}
\source{
TimiRGeN package.
}
\examples{
data(Clusters)
## maybe str(Clusters) ; plot(Clusters) ...
}
\keyword{datasets}
