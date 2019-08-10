library(data.table)
library(packHV) # required for plotting the histogram and boxplot together

file <- fread('gatk.vcf')
print (file)
quality <- file$QUAL
position <- file$POS

quality_histogram <- function(x){
  hist_boxplot(x,
       main = "Quality Distribution",
       xlab = "QUAL",
       col = 'cyan')
}
quality_histogram(quality)