#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(tidyverse)
#load data
dloadTargetScans() -> TargetScans
#checks
#check aspects of TargetScans
test_that("Check if this is a database of miR-mRNA predictions", {
expect_gt(length(unique(factor(TargetScans$ID))),
length(unique(factor(TargetScans$Transcript))))
expect_true(is.numeric(TargetScans$Symbol))
expect_equal(length(names(TargetScans)), 5)
})
#remove file
file.remove("Predicted_Targets_Context_Scores.default_predictions.txt")
#run function
TargetScans_data(targetScan = TargetScans, species = 'mmu') -> TargetScans_mmu
TargetScans_mmu[1:100,] -> TargetScans_mmu
#check 1
#should have three columns and be smaller than origional file
test_that("TargetScans_mmu is expected to be 3 columns long", {
expect_equal(length(names(TargetScans_mmu)), 3)
expect_gt(length(rownames(TargetScans)),
length(rownames(TargetScans_mmu)))
})
#internal checks
TargetScans %>%
filter(str_detect(Transcript, 'mmu')) -> TargetScans_s
#check 2
#TargetScans_s should be smaller than TargetScans
test_that("TargetScans_s is smaller than TargetScans", {
expect_equal(length(names(TargetScans_s)), 5)
expect_gt(length(rownames(TargetScans)),
length(rownames(TargetScans_s)))
})
#continue
tolower(TargetScans_s$ID) -> TargetScans_s$ID2
firstup <- function(x) {
substr(x, 1, 1) <- toupper(substr(x, 1, 1))
x
}
firstup(TargetScans_s$ID2) -> TargetScans_s$ID2
Targetscans_df <- data.frame(Targetscans_Interactions = paste(
TargetScans_s$Transcript,
':', TargetScans_s$ID2,sep = ''),
Targetscans_miR = TargetScans_s$Transcript,
Targetscans_mRNA = TargetScans_s$ID2)
Targetscans_df[1:100,] -> Targetscans_df
#check 2
#should be equal to function output
test_that("functional and manual output are the same", {
expect_equal(Targetscans_df, TargetScans_mmu)
})
#check 3
#test for human data
TargetScans_data(targetScan = TargetScans, species = 'hsa') -> TargetScans_hsa
#human data should be larger than mouse data
test_that("there are more human miR-mRNA predicted interactions than in
mice", {
expect_gt(length(rownames(TargetScans_hsa)),
length(rownames(TargetScans_mmu)))
})
#save data
as.matrix(Targetscans_df) -> TargetScans_matrix
saveRDS(TargetScans_matrix, "TargetScans_results.rds", compress = "xz")
