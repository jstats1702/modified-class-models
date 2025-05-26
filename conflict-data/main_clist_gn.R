# ------------------------------------------------------------
# Setup: Working directory, clean environment, load training function
# ------------------------------------------------------------

# Set working directory
setwd("/Users/juansosa/Dropbox/PAPERS/projects/eleni/conflict/")

# Clear environment
rm(list = ls())

# Load training function
source("training_function_clist_gn.R")

run_training(
     route   = "/Users/juansosa/Dropbox/PAPERS/projects/eleni/conflict/",
     data    = "conflict",
     model   = "clist", 
     prior   = "gn",
     n_burn  = 10000, 
     n_sams  = 100000, 
     n_skip  = 1
)