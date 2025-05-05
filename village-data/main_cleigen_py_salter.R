# ------------------------------------------------------------
# Setup: Working directory, clean environment, load training function
# ------------------------------------------------------------

# Set working directory
setwd("/Users/juansosa/Dropbox/PAPERS/projects/eleni/salter/")

# Clear environment
rm(list = ls())

# Load training function
source("training_function_cleigen_py.R")

run_training(
     route   = "/Users/juansosa/Dropbox/PAPERS/projects/eleni/salter/",
     data    = "salter",
     model   = "cleigen", 
     prior   = "py",
     n_burn  = 100000, 
     n_sams  = 1000000, 
     n_skip  = 1
)