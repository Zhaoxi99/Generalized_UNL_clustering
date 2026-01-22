rm(list=ls()) 
library(WASABI)
load("//csce.datastore.ed.ac.uk/csce/maths/groups/mdt/clustering_paper/application_data/dde/dde_fit_partitions.RData")
ptm=proc.time()
set.seed(129)
out_elbow <- elbow(fit_lddp$z, L_max = 6, psm = psm_lddp, ncores = parallel::detectCores()-2,
                   multi.start = 4, mini.batch = 300,
                   method.init = "++", method = "salso")
proc.time()-ptm
save(out_elbow,"//csce.datastore.ed.ac.uk/csce/maths/groups/mdt/clustering_paper/application_data/dde/dde_wasabi_elbow.RData")
