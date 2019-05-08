# Locations
file_name <- "z1_eec_var_global_80.csv"
out_folder <- "/Users/erinconrad/Desktop/residency stuff/R25/network subsets/results/for_r/"

# Add dunn test to library
library(multcomp)

# Load data
list_data <- read.csv(file=paste(out_folder,file_name,sep=""), header=FALSE, sep=",")
data <- matrix(unlist(list_data[,c("V1","V2","V3")]),nc=3)

# Do Friedman test
friedman.test(data)