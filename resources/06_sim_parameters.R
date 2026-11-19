library("dplyr")

gt_wf_params <- expand.grid(
    MU=7e-09,
    R=8.06452e-10,
    N=c(1000),
    G=10e6, # high enough so that more that there are at least L mutations post-burn-in
    L=c(1, 10, 100, 200),
    VA=c(0, 0.001, 0.01, 0.1),
    SIGMA=c(0, 0.1, 0.5, 0.9, 0.99),
    TAU=10,
    rep=1:30,
    adjust=c(T, F),
    # USE_VCF = 0 simulates a fresh genome and burns in de novo mutations
    # (gt_expectations.slim); USE_VCF = 1 seeds the run from standing variation in a
    # VCF window (gt_expectations_structured.slim). The USE_VCF = 1 rows are not
    # wired up yet -- that script still takes tmu/alpha -- so their shell branch
    # exits with an error; see workflow/rules/slim.smk.
    USE_VCF=0,
    type="wf"
)

gt_nonwf_params <- expand.grid(
    N = 1000,
    MU=7e-09,
    R=8.06452e-10,
    G=10e6, # simulated genome size when USE_VCF = 0
    SIGMA=c(0, 0.1, 0.5, 0.9, 0.99),
    VA=c(0, 0.001, 0.01, 0.1),
    L=c(1, 10, 100, 200),
    rep=1:30,
    TAU=10,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.25, 0.5, 0.75),
    BANK_SURV=c(0, 0.5, 0.75),
    type="nonwf",
    SURVIVAL_SELECTION=c(0, 0.5, 1),
    MIN_AGE=1,
    MAX_AGE=3,
    adjust = c(T,F),
    USE_VCF = 0
)

# VA = 0 makes alpha = 0, so neither L (the number of selected loci) nor
# SURVIVAL_SELECTION (how their variance is split above / below ground) has any
# effect: every combination yields the same neutral run. Keep one row per remaining
# combination and set both to NA rather than simulating the same design repeatedly.
#
# Only the nonWF blocks are collapsed. The wf grid is keyed by ID in the notebook
# (params[params$ID == simid, ]), and dropping rows there would renumber the 4,800
# wf simulations that already exist.
gt_nonwf_params <- gt_nonwf_params[gt_nonwf_params$VA > 0 |
                                     (gt_nonwf_params$L == min(gt_nonwf_params$L) &
                                        gt_nonwf_params$SURVIVAL_SELECTION == 0.5), ]
gt_nonwf_params$L[gt_nonwf_params$VA == 0] <- NA
gt_nonwf_params$SURVIVAL_SELECTION[gt_nonwf_params$VA == 0] <- NA

# adjust Ne based on selfing rate
gt_wf_params$N[(gt_wf_params$adjust == T)] = as.integer(gt_wf_params$N[(gt_wf_params$adjust == T)]/(1 - gt_wf_params$SIGMA[(gt_wf_params$adjust == T)]/2))

gt_nonwf_params$K[(gt_nonwf_params$adjust == T)] = as.integer(gt_nonwf_params$K[(gt_nonwf_params$adjust == T)]/(1 - gt_nonwf_params$SIGMA[(gt_nonwf_params$adjust == T)]/2))

gt_nonwf_params$BURN = 10*gt_nonwf_params$K

# combine sims
gt_params <- bind_rows(gt_wf_params, gt_nonwf_params)

# add simulation id
gt_params$ID <- 1:nrow(gt_params)

# save
write.table(gt_params, "../config/gt_params.tsv", sep = "\t", quote = F, row.names = F)
