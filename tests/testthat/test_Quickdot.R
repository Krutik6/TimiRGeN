#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(ggplot2)
#load data
readRDS("EnrichWiki.rds") -> sigwiki
#test Quickbar
Quickdot(X = sigwiki$D7_wikipathways,
         Y = sigwiki[4])

#test plot
ggplot()+
geom_dotplot(mapping = aes(x= reorder(Description, -p.adjust), y=Count,
fill=-p.adjust),
data = head(sigwiki$D7_wikipathways[which(
sigwiki$D7_wikipathways$p.adjust < 0.05),
], n = 15),
binaxis = 'y', dotsize = 2.5,
method = 'dotdensity', binpositions = 'bygroup', binwidth = 0.3,
stackdir = "center")+
scale_fill_continuous(type = "gradient") +
labs(y = "Associated genes", x = "wikipathways", fill = "p.adjust") +
theme(axis.text=element_text(size=8)) +
ggtitle(names(sigwiki[4])) +
theme(plot.title = element_text(2, face = "bold", hjust = 0.5),
legend.key.size = unit(2, "line")) +
theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
coord_flip()
#visual check
test_that("should be 4 wikipathways", {
  expect_equal(length(rownames(sigwiki$D7_wikipathways@result)),4)})
#outputs are the same
