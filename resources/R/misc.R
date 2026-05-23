box::use(dplyr[...])
box::use(utils[write.table])

#' Little simulation, testing if no-multicollinearity == unbiased estimates
#'
#' @param n 
#' @param a 
#' @param b 
#'
#' @return
#' @export
#'
#' @examples hist(replicate(1000, lm_sim(100, -10, 10)))
lm_sim <- function(n, a, b) {
  
  x1 <- sample(1:100, size = n, replace = T)
  x2 <- sample(1:100, size = n, replace = T)
  
  e <- rnorm(n)
  
  y = a*x1 + b*x2 + e
  
  coefficients(lm(y ~ x1))[2]
}


#' Simulation to test how multi-collinearity affects estimates
#'
#' @param mu1 
#' @param mu2 
#' @param covmat 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
multicol_sim <- function(n, mu1, mu2, b0, b1, b2, covmat){
  
  # Generate 1,000 random samples
  bivariate_data <- mvrnorm(n = n, mu = c(mu1, mu2), Sigma = covmat)
  
  # Convert to a data frame for easier use
  df <- as.data.frame(bivariate_data)
  colnames(df) <- c("X1", "X2")
  
  # make linear model
  epsilon <- rnorm(n)
  df$Y <- b0 + b1*df$X1 + b2*df$X2 + epsilon
  
  # estimate parameters
  modf <- melt(coefficients(lm(Y ~ X2 + X1, data = df)))
  modf$type <- "joint"
  
  mod1 <- melt(coefficients(lm(Y ~ X1, data = df)))
  mod1$type <- "single"
  
  mod2 <- melt(coefficients(lm(Y ~ X2, data = df)))
  mod2$type <- "single"
  
  all_mod <- rbind(modf, mod1, mod2)
  all_mod$term <- gsub("X11", "X1", rownames(all_mod))
  all_mod$term <- gsub("X21", "X2", all_mod$term)
  all_mod$term <- gsub(").", ")", all_mod$term)
  
  return(all_mod)
}


#' Mark entries after first NA as also NA
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
rm_na_after_na <- function(x){
  
  # get first na in sequence
  first_na <- min(which(is.na(x)))
  
  # mark instances after first na as also na
  x[first_na:length(x)] <- NA
} 


#' Append row to table, good for for loops
#'
#' @param x 
#' @param output_name 
#'
#' @return
#' @export
#'
#' @examples
append_table <- function(x, output_name) {
  write.table(x, 
              output_name, 
              col.names=!file.exists(output_name), 
              append = T, 
              row.names = F, 
              sep = ",", 
              quote = F)
}

#' Block data by chromosome and position
#'
#' @param data 
#' @param chrom_col_name 
#' @param pos_col_name 
#' @param window_size 
#' @param sep 
#'
#' @return
#' @export
#'
create_blocks <- function(data, chrom_col_name, pos_col_name, window_size, sep){

  sort_check <- data %>%
    group_by({{ chrom_col_name }}) %>%
    summarise(is_sorted = all(diff({{ pos_col_name }}) >= 0))

  if(!all(sort_check$is_sorted)){
    stop("Input data should be sorted by position within chromosome.")
  }

  blocks <- data %>%
    group_by({{ chrom_col_name }}) %>%
    mutate(block = (row_number() - 1) %/% window_size + 1,
           window = paste({{ chrom_col_name }}, block, sep = sep))

  return(blocks)
}

#' Grab sample sizes for a subset of samples
#'
#' @param n_data 
#' @param sample_col 
#' @param samples 
#' @param first_n
#'
#' @return
#' @export
#'
grab_sample_sizes <- function(n_data, sample_col, samples, first_n){
    rep_n <- n_data %>% filter({{ sample_col }} %in% samples)
    rep_n <- rbind(c(0, first_n), rep_n)
    rep_n <- rep_n %>% arrange({{ sample_col }})
    return(rep_n)
}

#' Fit allele frequency under selection
#'
#' @param q 
#' @param h 
#' @param s 
#'
#' @return
#' @export
#'
fitfreq <- function(q, h, s){
  p=1-q;
  return((q^2*(1+s) + p*q*(1+h*s))/( 1 + s*q*(2*h*p+q)))
}

WF.sel <- function(N, q, h, s, G){
  t=array(,dim=G)
  t[1] = N*q
  for(i in 2:G){
    t[i] = rbinom(1,N,fitfreq(t[i-1]/N, h, s[i-1]))
  }
  return(t)
}

check_dp_match_dt <- function(x,y){
  results = NULL
  for(i in 1:nrow(x)){
    mean_cov <- x[i,"mean_cov"]
    int1 <- x[i,"Var1"]
    int2 <- x[i,"Var2"]
    site <- as.numeric(x[i,"site"])

    dT1 <- y[(y$site == site) & (y$interval == int1),"dT"]

    dT2 <- y[(y$site == site) & (y$interval == int2),"dT"]

    if(dT1*dT2 > 0){
      results <- rbind(results, data.frame(site = site, int1 = int1, int2 = int2, match = T))
    } else {
      results <- rbind(results, data.frame(site = site, int1 = int1, int2 = int2, match = F))
    }
  }
  print(paste("Done at : ", Sys.time()))
  return(table(results$match))
}

