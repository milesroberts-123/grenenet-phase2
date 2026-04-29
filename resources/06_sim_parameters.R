library("dplyr")

gt_wf_params <- expand.grid(
    nmu=7e-09,
    tmu=c(7e-10),
    R=8.06452e-10,
    N=c(1000),
    L=1e6,
    sigma=c(0, 0.05, 0.5, 0.95, 0.99),
    alpha=c(0, 0.0025, 0.005, 0.01, 0.015),
    gamma=5,
    tau=105,
    rep=1:30,
    adjust=c(T),
    #type="bank",
    type=c("unstruct", "struct")
)

gt_nonwf_params <- expand.grid(N = 1000,
			       nmu=7e-09,
			       tmu=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    sigma=c(0, 0.05, 0.5, 0.95, 0.99),
    alpha=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    gamma=10,
    tau=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.9),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    MIN_AGE=1,
    MAX_AGE=20,
    adjust = T
)

gt_nonwf_wflike <- expand.grid(N = 1000,
                               nmu=7e-09,
                               tmu=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    sigma=c(0, 0.05, 0.5, 0.95, 0.99),
    alpha=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    gamma=10,
    tau=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.0),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    adjust = T,
    MIN_AGE=0,
    MAX_AGE=20 # doesn't mattern because BANK_SURV=0
			       ) 

gt_nonwf_params <- rbind(gt_nonwf_params, gt_nonwf_wflike)

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$sigma[(gt_wf_params$adjust == T)]/2))

# combine sims

# add simulation id
gt_nonwf_params$ID <- 1:nrow(gt_nonwf_params)

# save
write.table(gt_nonwf_params, "../config/gt_params.tsv", sep = "\t", quote = F, row.names = F)

# fst simulations parameters
fst_params <- expand.grid(
    mu=7e-09,
    R=8.06452e-10,
    msprimeN=24000,
    L=c(3081,6583),
    sigma=0.95,
    rep=1:10000,
    tau=14)

# adjust Ne based on selfing rate
fst_params$slimN = as.integer(fst_params$msprimeN/(1 - fst_params$sigma/2))

# add simulation id
fst_params$ID <- 1:nrow(fst_params)

# save
write.table(fst_params, "../config/fst_params.tsv", sep = "\t", quote = F, row.names = F)
