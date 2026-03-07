library("dplyr")

params <- expand.grid(
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
params$N[(params$adjust == T)] = as.integer(params$N[(params$adjust == T)]/(1 - params$sigma[(params$adjust == T)]/2))

# add simulation id
params$ID <- 1:nrow(params)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
