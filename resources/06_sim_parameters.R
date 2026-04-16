library("dplyr")

gt_params <- expand.grid(
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
    struct="struct"
    #struct=c("unstruct", "struct")
)

# adjust Ne based on selfing rate
gt_params$N[(gt_params$adjust == T)] = as.integer(gt_params$N[(gt_params$adjust == T)]/(1 - gt_params$sigma[(gt_params$adjust == T)]/2))

# add simulation id
gt_params$ID <- 1:nrow(gt_params)

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
