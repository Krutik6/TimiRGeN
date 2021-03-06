% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiReg.R
\name{multiReg}
\alias{multiReg}
\title{multiReg}
\usage{
multiReg(MAE, gene_interest, mRNAreg, filt_df, miRNA_exp, mRNA_exp)
}
\arguments{
\item{MAE}{Input MAE which stores results from multiReg It is
recommended to use the MAE which was used in matrixFilter.}

\item{gene_interest}{Name of gene of interest. Full string required with no
spelling mistakes e.g. "mmu-miR-140-5p", "Ffg1", "hsa-miR-29a-5p". If gene
of interest is an mRNA, mRNAreg must be TRUE. Otherwise, if the gene of interest
is a miRNA, then mRNAreg must be FALSE.}

\item{mRNAreg}{TRUE or FALSE. Is the gene of interest a mRNA? Default is TRUE.}

\item{filt_df}{Dataframe from the matrixFilter function.}

\item{miRNA_exp}{miRNA data from using the diffExpressRes function on miRNA
data.}

\item{mRNA_exp}{mRNA data from using the diffExpressRes function on miRNA
data}
}
\value{
A matrix which contains the gene of interest and all binding partners.
Their values (Log2FC or ave exp) for each time point are also produced.
}
\description{
Creates a new matrix containing the gene of interest and
each binding partner that it interacts with.
}
\examples{
library(org.Mm.eg.db)

miR <- mm_miR[1:100,]

mRNA <- mm_mRNA[1:200,]

MAE <- startObject(miR = miR, mRNA = mRNA)

MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')

MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus', )

MAE <- diffExpressRes(MAE, df = assay(MAE, 1), dataType = 'Log2FC',
                     genes_ID = assay(MAE, 3),
                      idColumn = 'GENENAME',
                      name = "miRNA_log2fc")

MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
                     genes_ID = assay(MAE, 7),
                     idColumn = 'GENENAME',
                     name = "mRNA_log2fc")

Filt_df <- data.frame(row.names = c("mmu-miR-145a-3p:Adamts15",
                                   "mmu-miR-146a-5p:Acy1"),
                     corr = c(-0.9191653, 0.7826041),
                     miR = c("mmu-miR-145a-3p", "mmu-miR-146a-5p"),
                     mRNA = c("Adamts15", "Acy1"),
                     miR_Entrez = c(387163, NA),
                     mRNA_Entrez = c(235130, 109652),
                     TargetScan = c(1, 0),
                     miRDB = c(0, 0),
                     Predicted_Interactions = c(1, 0),
                     miRTarBase = c(0, 1),
                     Pred_Fun = c(1, 1))

MAE <- matrixFilter(MAE, miningMatrix = Filt_df, negativeOnly = FALSE,
                   threshold = 1, predictedOnly = FALSE)

MAE <- multiReg(MAE = MAE, gene_interest = "Adamts15",
                mRNAreg =TRUE, filt_df=MAE[[11]], miRNA_exp=MAE[[9]],
                mRNA_exp=MAE[[10]])
}
