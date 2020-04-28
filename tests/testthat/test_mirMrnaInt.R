library(clusterProfiler)
library(org.Mm.eg.db)
library(testthat)
#load data
#test function
G <- data.frame(row.names = c("Acaa1a", "Acadm", "Acss1", "Adh1"),
                "D1.Log2FC" = c("-1.2944593","-2.0267432","-2.1934942",
                                "-2.1095853"),
                "D2.Log2FC" = c("-1.1962396","-2.1345451","-1.7699232",
                                "-1.0961674"),
                "D3.Log2FC" = c("0.2738496","-1.9991046","-1.7637549",
                                "-1.6572653"),
                "D7.Log2FC" = c("-0.51765245","-2.20689661","-0.68479699",
                                "-2.06512466"),
                "D14.Log2FC" = c("-0.4510294","-1.1523849","-0.4297012",
                                 "-1.1017597"),
                "ID" = c("113868","11364","68738","11522"))

MIR <- data.frame(row.names = c("mmu-miR-101a-3p", "mmu-miR-101a-5p",
                                "mmu-miR-101c", "mmu-miR-106a-5p"),
                  "D1.Log2FC" = c("-0.0039141722","-0.4328659746",
                                  "-0.0038897133", "-0.4161749123"),
                  "D2.Log2FC" = c("-0.210605345","-0.600422732",
                                  "-0.210574742", "-0.530311376"),
                  "D3.Log2FC" = c("-0.315070839","-0.745367163",
                                  "-0.315012148", "-0.559274530"),
                  "D5.Log2FC" = c("-0.41087763","-0.63952382",
                                  "-0.41087876", "-1.03618015"),
                  "D14.Log2FC" = c("-0.39466968","-0.60122678",
                                   "-0.39461099", "-0.41889698"),
                  "ID" = c("387143","387143","100628572","723829"))

MAE <- MultiAssayExperiment(list(G = G, MIR = MIR))

allnames <- expand.grid(rownames(assay(MAE, 2)), rownames(assay(MAE, 1)),
                        stringsAsFactors = FALSE)

interactions_df <- do.call(rbind, apply(allnames,
                                        1,
                                        corrTable,
                                        assay(MAE, 2),
                                        assay(MAE, 1),
                                        maxInt = 5))

test_that("expect 5 rows and 16 columns (4*4)", {
    expect_equal(length(rownames(interactions_df)), 16)
    expect_equal(length(colnames(interactions_df)), 5)
})

MAE2 <- suppressMessages(MultiAssayExperiment(list(
                         interactions_df = interactions_df)))
MAE <- c(MAE, MAE2)
