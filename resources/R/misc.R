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

grab_sample_sizes <- function(n_data, sample_col, samples, first_n){
    rep_n <- n_data %>% filter({{ sample_col }} %in% samples)
    rep_n <- rbind(c(0, first_n), rep_n)
    rep_n <- rep_n %>% arrange({{ sample_col }})
    return(rep_n)
}

fitfreq = function(q, h, s){
  p=1-q;
  return((q^2*(1+s) + p*q*(1+h*s))/( 1 + s*q*(2*h*p+q)))
}

WF.sel=function(N, q, h, s, G){
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

