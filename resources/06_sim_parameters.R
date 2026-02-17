library("dplyr")

params <- expand.grid(
    nmu=1e-8,
    tmu=2e-9,
    R=1e-8,
    N=1000,
    L=5e6,
    n=50,
    sigma=c(0, 0.1, 0.5, 0.9),
    alpha=c(0, 0.005, 0.01),
    gamma=5,
    tau=95,
    rep=1:30
)

# add simulation id
params$ID <- 1:nrow(params)

# save
write.table(params, "../config/parameters.tsv", sep = "\t", quote = F, row.names = F)
