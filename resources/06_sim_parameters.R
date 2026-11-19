library("dplyr")

gt_wf_params <- expand.grid(
    MU=7e-09,
    R=8.06452e-10,
    N=c(1000),
    G=10e6, # high enough so that more that there are at least L mutations post-burn-in
    L=c(1, 10, 100, 200),
    VA=c(0, 0.001, 0.01, 0.1),
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    TAU=10,
    rep=1:30,
    adjust=c(T, F),
    # USE_VCF = 0 simulates a fresh genome and burns in de novo mutations
    # (gt_expectations.slim); USE_VCF = 1 seeds the run from standing variation in a
    # VCF window (gt_expectations_structured.slim). The USE_VCF = 1 rows are not
    # wired up yet -- that script still takes tmu/alpha -- so their shell branch
    # exits with an error; see workflow/rules/slim.smk.
    USE_VCF=c(0, 1),
    type="wf"
)

gt_nonwf_bank_struc <- expand.grid(
    N = 1000,
    MU=7e-09,
    R=8.06452e-10,
    G=10e6, # simulated genome size when USE_VCF = 0
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    VA=c(0, 0.001, 0.01, 0.1),
    L=c(1, 10, 100, 200),
    rep=1:30,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.15, 1),
    BANK_SURV=c(0.5, 0.8, 1),
    type="bank",
    SURVIVAL_SELECTION=c(0.5,1),
    MIN_AGE=1,
    MAX_AGE=3,
    adjust = T,
    USE_VCF = 0
)

gt_nonwf_bank_nostruc <- expand.grid(N = 1000,
                               MU=7e-09,
                               R=8.06452e-10,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    VA=c(0, 0.001, 0.01, 0.1),
    G=10e6, # simulated genome size when USE_VCF = 0
    L=c(1, 10, 100, 200),
    rep=1:30,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.6),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    MIN_AGE=1,
    MAX_AGE=3,
    adjust = T,
    USE_VCF = 0
)

gt_nonwf_nobank_struc <- expand.grid(N = 1000,
                               MU=7e-09,
                               R=8.06452e-10,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    VA=c(0, 0.001, 0.01, 0.1),
    G=10e6, # simulated genome size when USE_VCF = 0
    L=c(1, 10, 100, 200),
    rep=1:30,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.0),
    type="bank",
    SURVIVAL_SELECTION=c(0.5, 1),
    adjust = T,
    MIN_AGE=0, # no seed bank (BANK_SURV = 0): seeds must germinate immediately
    MAX_AGE=3, # irrelevant while BANK_SURV = 0
    USE_VCF=0
    ) 

#gt_nonwf_params <- rbind(gt_nonwf_params, gt_nonwf_wflike)

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$SIGMA[(gt_wf_params$adjust == T)]/2))

gt_nonwf_nobank_nostruc <- expand.grid(N = 1000,
                               MU=7e-09,
                               R=8.06452e-10,
    SIGMA=c(0, 0.05, 0.5, 0.95, 0.99),
    VA=c(0, 0.001, 0.01, 0.1),
    G=10e6, # simulated genome size when USE_VCF = 0
    L=c(1, 10, 100, 200),
    rep=1:30,
    TAU=12,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.8),
    BANK_SURV=c(0.0),
    type="bank",
    SURVIVAL_SELECTION=0.5,
    adjust = T,
    MIN_AGE=0,
    MAX_AGE=3, # doesn't matter because BANK_SURV=0
    USE_VCF=0
)

gt_nonwf_params <- bind_rows(gt_nonwf_bank_struc, gt_nonwf_bank_nostruc, gt_nonwf_nobank_struc, gt_nonwf_nobank_nostruc)

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$SIGMA[(gt_wf_params$adjust == T)]/2))

gt_nonwf_params$K[(gt_nonwf_params$adjust == T)] = as.integer(gt_nonwf_params$K[(gt_nonwf_params$adjust == T)]/(1 - gt_nonwf_params$SIGMA[(gt_nonwf_params$adjust == T)]/2))

gt_nonwf_params$BURN = 10*gt_nonwf_params$K

# combine sims
gt_params <- bind_rows(gt_wf_params, gt_nonwf_params)

# add simulation id
gt_params$ID <- 1:nrow(gt_params)

# which simulation families enter config/gt_params.tsv
sim_types <- c("wf", "bank")
gt_params <- gt_params %>% filter(type %in% sim_types)

# save
write.table(gt_params, "../config/gt_params.tsv", sep = "\t", quote = F, row.names = F)
