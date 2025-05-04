setwd("C:/Users/User/Dropbox/PAPERS/projects/eleni")

load("simulation_class_dm_hom_equal_medium.RData")
plot (samples$loglik_chain, type = "p", pch = ".", col = adjustcolor(1, 0.3))

load("simulation_class_dp_hom_equal_medium.RData")
lines(samples$loglik_chain, type = "p", pch = ".", col = adjustcolor(2, 0.3))

load("simulation_class_py_hom_equal_medium.RData")
lines(samples$loglik_chain, type = "p", pch = ".", col = adjustcolor(3, 0.3))
