library("dplyr")

gt_wf_params <- expand.grid(
    MU=7e-09,
    R=8.06452e-10,
    N=c(1000),
    G=5e6, # high enough so that more that there are at least L mutations post-burn-in
    L=c(10, 100, 1000),
    VA=c(0.01, 0.05, 0.1),
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    TAU=10,
    rep=1:30,
    adjust=c(T, F),
    type=c("unstruct", "struct")
)

gt_nonwf_bank_struc <- expand.grid(
    N = 1000,
    MU=7e-09,
    TMU=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    ALPHA=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    GAMMA=10,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.15, 1),
    BANK_SURV=c(0.5, 0.8, 1),
    type="bank",
    SURVIVAL_SELECTION=c(0.5,1),
    MIN_AGE=1,
    MAX_AGE=5,
    adjust = T,
    STRUCTURE = 1,
    BURN=999
)

gt_nonwf_bank_nostruc <- expand.grid(N = 1000,
                               MU=7e-09,
                               TMU=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    ALPHA=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    GAMMA=10,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.6),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    MIN_AGE=1,
    MAX_AGE=5,
    adjust = T,
    STRUCTURE = 0,
    BURN=10000
)

gt_nonwf_nobank_struc <- expand.grid(N = 1000,
                               MU=7e-09,
                               TMU=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    ALPHA=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    GAMMA=10,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.0),
    type="bank",
    SURVIVAL_SELECTION=c(0.5, 1),
    adjust = T,
    MAX_AGE=20 # doesn't mattern because BANK_SURV= 0
    ) 

#gt_nonwf_params <- rbind(gt_nonwf_params, gt_nonwf_wflike)

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$SIGMA[(gt_wf_params$adjust == T)]/2))

gt_nonwf_nobank_nostruc <- expand.grid(N = 1000,
                               MU=7e-09,
                               TMU=c(2e-11),
    R=8.06452e-10,
    L=1e6,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    ALPHA=c(0, 0.0025, 0.005, 0.01, 0.015),
    rep=1:30,
    GAMMA=10,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.0),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    adjust = T,
    MIN_AGE=0,
    MAX_AGE=5, # doesn't mattern because BANK_SURV=0
    STRUCTURE=0,
    BURN=10000    )

gt_nonwf_params <- bind_rows(gt_nonwf_bank_struc, gt_nonwf_bank_nostruc, gt_nonwf_nobank_struc, gt_nonwf_nobank_nostruc)

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$SIGMA[(gt_wf_params$adjust == T)]/2))

gt_nonwf_params$K[(gt_nonwf_params$adjust == T)] = as.integer(gt_nonwf_params$K[(gt_nonwf_params$adjust == T)]/(1 - gt_nonwf_params$SIGMA[(gt_nonwf_params$adjust == T)]/2))

gt_nonwf_params$BURN = 10*gt_nonwf_params$K

# combine sims
gt_params <- bind_rows(gt_wf_params, gt_nonwf_params)

# add simulation id
gt_params$ID <- 1:nrow(gt_params)

# subset only 
gt_params <- gt_params %>% filter(type == "unstruct")

# save
write.table(gt_params, "../config/gt_params.tsv", sep = "\t", quote = F, row.names = F)

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
