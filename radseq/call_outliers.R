# Data file and thresholds are provided as a command line argument to the script.
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
residuals_sd_cutoff <- -as.numeric(args[2])

# Validate arguments.
if (!is.character(filename) || !file.exists(filename)) quit(status(1))
if (!is.numeric(residuals_sd_cutoff)) quit(status(1))

# Read table and add residuals from regressing percent homozygous calls to
# percent called as the fourth column. 
data <- read.table(filename, sep = "\t", header = FALSE)
data$V4 = rstandard(lm(V3 ~ V2 + 0, data = data))

# Outliers are the individuals with very few called sites or high residual.
outliers = subset(data, V4 < residuals_sd_cutoff)
for (name in outliers[,1]) {
    cat(name, sep ="\n");
}