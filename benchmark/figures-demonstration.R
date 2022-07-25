# The following experiments compare different algorithms for approximate inference

# Load packages

library("SGS")
library("usethis")

# packages for benchmark visualization

library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("cowplot")
library("gridExtra")
library("matrixStats")

###############################
## Test Different Dimensions ##
###############################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkDims.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

p1 <- SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("50","100","150","200"), 
             c("Nodes N=50","Nodes N=100","Nodes N=150","Nodes N=200"), "Dims2", width = 5, height = 4, niceBoxPlot=TRUE, retP=TRUE)

#########################
## Test Different Nets ##
#########################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkNets.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

# final plots 

p2 <- SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("E-R","E-R Isl.","B-A","W-S"), 
             c("Erdoes-Renyi graph","Erdoes-Renyi Island graph","Barabási-Albert graph","Watts-Strogatz graph"), "Nets2", width = 5, height = 4, niceBoxPlot=TRUE, retP=TRUE)

#############################
## Test Different Sparcity ##
#############################

# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkSparcity.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

p3 <- SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("S=2","S=3","S=4","S=5"), 
             c("Markov Blanket size S=2","Markov Blanket size S=3","Markov Blanket size S=4","Markov Blanket size S=5"), "Sparcity2", width = 5, height = 4, niceBoxPlot=TRUE, retP=TRUE)

#####################################
## Different Fractions of Evidence ##
#####################################

# # load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkFractions.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

p4 <- SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("f=0.2","f=0.4","f=0.6","f=0.8"), 
             c("Evidence Fraction f=0.2","Evidence Fraction f=0.4","Evidence Fraction f=0.6","Evidence Fraction f=0.8"), "Fractions2", width = 5, height = 4, niceBoxPlot=TRUE, retP=TRUE)

#########################
## Multiple Categories ##
#########################

# load benchmark results
benchmarkResutls <- readRDS(file = "benchmark_res/benchmarkCategories.rds")
resultTable1 <- benchmarkResutls$resultTable1
resultTable2 <- benchmarkResutls$resultTable2
resultTable3 <- benchmarkResutls$resultTable3
resultTable4 <- benchmarkResutls$resultTable4

## FinalPlots:

p5 <- SGS:::makeAllPlots(resultTable1,resultTable2,resultTable3, resultTable4, c("C=2","C=4","C=6", "C=8"), 
             c("Categories C=2","Categories C=4","Categories C=6","Categories C=8"), "Categories2", width = 5, height = 4, niceBoxPlot=TRUE, retP=TRUE)


###################
## Summary Plots ##
###################

pall <- cowplot::plot_grid(p1+theme(legend.position="none")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                   p2+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.y=element_blank(),legend.position="none")+labs(x="Graph type"),
                   p3+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.y=element_blank(),legend.position="none")+labs(x="Markov blanket size"),
                   p4+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.y=element_blank(),legend.position="none")+labs(x="Evidence fraction"),
                   p5+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.y=element_blank(),legend.position="none")+labs(x="Category size"), 
                   ncol = 5,rel_widths = c(1.14,1,1,1,1))
legend <- cowplot::get_legend(p1)

p6s <- cowplot::plot_grid(pall, legend,nrow = 2, rel_heights = c(1, 0.1), rel_widths = c(1,1))

# cairo_pdf(paste0("~/Desktop/summary.pdf"), width = 8.5, height = 3.7)
# p6s
# dev.off()

