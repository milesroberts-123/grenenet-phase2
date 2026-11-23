library("dplyr")

# The wf grid, for a given set of selfing rates. Used twice below so that the
# original selfing levels keep their row positions: the notebook joins parameters
# to results by ID (params[params$ID == simid, ]), so the wf simulations that
# already exist must not move.
gt_wf_params <- expand.grid(
        MU=7e-09,
        R=8.06452e-10,
        N=c(1000),
        G=5e6, # genome size: simulated genome (wf), founder VCF window (nonWF)
        L=c(1, 10, 100),
        VA=c(0, 0.001, 0.01, 0.1),
        SIGMA=c(0, 0.1, 0.5, 0.95),
        TAU=10,
        rep=1:30,
        adjust=c(T, F),
        USE_VCF=0, # no VCF: the wf runs simulate a fresh genome and burn in
        type="wf"
    )

gt_nonwf_params <- expand.grid(
    N = 1000,
    MU=7e-09,
    R=8.06452e-10,
    G=5e6, # genome size: simulated genome (wf), founder VCF window (nonWF)
    SIGMA=c(0, 0.1, 0.5, 0.95),
    VA=c(0, 0.001, 0.01, 0.1),
    L=c(1, 10, 100),
    rep=1:30,
    TAU=10,
    K=c(1000),
    N_OFFSPRING=c(5),
    GERM_RATE=c(0.5, 0.75),
    BANK_SURV=c(0, 0.5),
    type="nonwf",
    SURVIVAL_SELECTION=c(0, 0.5, 1),
    MIN_AGE=1,
    MAX_AGE=3,
    adjust=T,
    USE_VCF=1
)

# BANK_SURV = 0 is the no-seed-bank control. MIN_AGE = 1 would make every seed die
# in storage before it becomes eligible to germinate (the population goes extinct
# within two ticks), so those rows germinate in their birth tick and run as a
# strict annual.
gt_nonwf_params$MIN_AGE[gt_nonwf_params$BANK_SURV == 0] <- 0

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
