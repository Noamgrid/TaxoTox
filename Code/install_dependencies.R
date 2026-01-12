# Install script for TaxoTox Shiny App dependencies
# This script checks if packages are already installed before trying to install them.

# List of required packages from CRAN
required_packages <- c(
    "shiny", 
    "openxlsx", 
    "tidyr", 
    "dplyr", 
    "ggplot2", 
    "scales", 
    "ggthemes", 
    "tidyverse", 
    "ggpattern", 
    "ggpubr", 
    "ggpmisc", 
    "hrbrthemes", 
    "remotes", 
    "purrr", 
    "readxl", 
    "stringr", 
    "RSQLite", 
    "data.table", 
    "stringdist", 
    "fst",
    "DT"
)

# Identify which packages are not yet installed
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, repos = "http://cran.us.r-project.org")
} else {
    cat("All required packages are already installed.\n")
}

cat("Dependency check complete.\n")

