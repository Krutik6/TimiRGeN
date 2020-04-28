#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load data
net <- readRDS("net.rds")
quickNet(metadata(net)[[1]])
dev.off()
#visual check, legends a bit wonky
test_that("net should be 10 long", {
  expect_equal(length(metadata(net)[[1]]), 10)
})
