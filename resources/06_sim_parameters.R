library("dplyr")

gt_params <- expand.grid(
    nmu=1e-8,
    tmu=1e-9,
    R=1e-8,
    N=1000,
    L=5e6,
    sigma=c(0, 0.05, 0.5, 0.95, 0.99),
    alpha=c(0, 0.01),
    gamma=5,
    tau=105,
    rep=1:60,
    adjust=c(T, F)
)

# adjust Ne based on selfing rate
gt_params$N[(gt_params$adjust == T)] = as.integer(gt_params$N[(gt_params$adjust == T)]/(1 - gt_params$sigma[(gt_params$adjust == T)]/2))

# add simulation id
gt_params$ID <- 1:nrow(gt_params)

# save
write.table(gt_params, "../config/gt_params.tsv", sep = "\t", quote = F, row.names = F)

# fst simulations parameters
fst_params <- expand.grid(
    mu=1e-8,
    R=1e-8,
    N=1000,
    L=5e6,
    sigma=c(0.95),
    rep=1:2000,
    adjust=c(T, F)
)

# adjust Ne based on selfing rate
fst_params$N[(fst_params$adjust == T)] = as.integer(fst_params$N[(fst_params$adjust == T)]/(1 - fst_params$sigma[(fst_params$adjust == T)]/2))

# add simulation id
fst_params$ID <- 1:nrow(fst_params)

# save
write.table(fst_params, "../config/fst_params.tsv", sep = "\t", quote = F, row.names = F)
