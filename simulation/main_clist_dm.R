# Set working directory (ensure correct path is used)
setwd("C:/Users/User/Dropbox/PAPERS/projects/eleni")  # setwd("C:/Users/Juan Camilo/Dropbox/PAPERS/projects/eleni")

source("simulation_function_clist_dm.R")

# --- Medium ---

run_simulation(model = "clist",
               prior = "dm",
               density_type = "hom",
               cluster_type = "equal",
               size = "medium",
               n_burn = 100000,
               n_sams = 1000000,
               n_skip = 1)

run_simulation(model = "clist",
               prior = "dm",
               density_type = "het",
               cluster_type = "random",
               size = "medium",
               n_burn = 100000,
               n_sams = 1000000,
               n_skip = 1)

# --- Large ---

run_simulation(model = "clist",
               prior = "dm",
               density_type = "het",
               cluster_type = "random",
               size = "large",
               n_burn = 100000,
               n_sams = 1000000,
               n_skip = 1)

# --- metrics ---