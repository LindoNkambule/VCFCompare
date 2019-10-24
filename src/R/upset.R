#!/usr/bin/env Rscript
#
# Plot VCFCompare.py UpSet Plots
#
# This script requires UpSetR and argepase packages.
# To install the packages, run the following commandc in R:
#
#   install.packages("UpSetR")
#   install.packages("argeparse")
#
# Usage:
#
# This script runs on a set of VCFCompare.py results. Result n will have been run
# with VCFCompare.py -o prefix_n.
#
# run Rscript upsetr.Rscript [-i] file.csv -o output_prefix

library(UpSetR)
library(argparse)

#Command line
parser <- ArgumentParser()
parser <- ArgumentParser(description='This is an extension of VCFCompare. Use this script to plot upset plots.')
parser$add_argument('-i', help='This is the input CSV file generated using VCFCompare.py')
parser$add_argument('-o', help='Output prefix filename for upset plots')
args <- parser$parse_args()
input <- args$i
output_prefix <- args$o

#Error handling
if (!file.exists(input)){
sprintf("File %s does not exist.", input)
stop(input, " does not exist")
}

file = read.csv(input)

### 1.SNPs ###
truthSNVunique <- file$FN[1]
querySNVunique <- file$FP[1]
SNVintersect <- file$TP[1]
expressionSNV <- c(Truth = truthSNVunique, Query = querySNVunique,
                        `Truth&Query` = SNVintersect)

pdf(paste(output_prefix, "SNV", "pdf", sep="."))
upset(fromExpression(expressionSNV),
      main.bar.color = "black",
      sets.bar.color = "black",
      matrix.color = "black",
      point.size=5,
      order.by = "freq")
dev.off()

### 2.INDELs ###
truthINDELunique <- file$FN[2]
queryINDELunique <- file$FP[2]
INDELintersect <- file$TP[2]
expressionINDEL <- c(Truth = truthINDELunique, Query = queryINDELunique,
                   `Truth&Query` = INDELintersect)

pdf(paste(output_prefix, "INDEL", "pdf", sep="."))
upset(fromExpression(expressionINDEL),
      main.bar.color = "black",
      sets.bar.color = "black",
      matrix.color = "black",
      point.size=5,
      order.by = "freq")
dev.off()
