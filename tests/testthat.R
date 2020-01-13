#unit test check
setwd("~/Documents/Package/smiRk/tests/testthat/")
library(testthat)
library(smiRk)
list.files() -> X
for (i in 1:length(X)) {
test_dir(X[i])
}
#bioconductor check
#library(BiocCheck)
#BiocCheck("~/Documents/Package/smiRk/")
