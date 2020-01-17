#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
library(org.Mm.eg.db)
library(clusterProfiler)
mm_miR -> miR
getIDs_miR_mouse(miR)
Filt_df <- data.frame(row.names = c("mmu-miR-320-3p:Acss1",
"mmu-miR-27a-3p:Odc1"),
avecor = c(-0.9191653, 0.7826041),
miR = c("mmu-miR-320-3p", "mmu-miR-27a-3p"),
mRNA = c("Acss1", "Acss1"),
miR_Entrez = c(NA, NA),
mRNA_Entrez = c(68738, 18263),
TargetScan = c(1, 0),
miRDB = c(0, 0),
Predicted_Interactions = c(1, 0),
miRTarBase = c(0, 1),
Pred_Fun = c(1, 1))
MakeMapp(filt_df = Filt_df, miR_IDs_adj = miR_adjusted_entrez,
Datatype = 'L') -> MAPPdata
#check 1
#check aspects of MAPPdata
test_that("MAPPdata runs as expected", {
expect_equal(length(rownames(Filt_df)),
length(rownames(MAPPdata)))
expect_equal(length(names(MAPPdata)), 3)
})
#check 2
#check if 'En' works
MakeMapp(filt_df = Filt_df, miR_IDs_adj = miR_adjusted_ensembl,
Datatype = 'En') -> MAPPdata2
test_that("En and L lead to different results", {
expect_equal(MAPPdata[,1], MAPPdata2[,1])
expect_false(isTRUE(all.equal(MAPPdata, MAPPdata2)))
})
#save
saveRDS(MAPPdata, "MAPPdata.rds", compress = "xz")

