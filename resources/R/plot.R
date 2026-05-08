box::use(ggplot2[...])
box::use(reshape2[melt])
box::use(dplyr[...])

get_upper_tri <- function(mat){
  mat[lower.tri(mat)]<- NA
  diag(mat) <- NA
  return(mat)
}

plot_var_cov_matrix <- function(covmat, output_name, include_values = TRUE, low_col = "#65014B", mid_col = "#F5F0F0", high_col = "#0C4C00", legend_name = "Cov(\u0394pi, \u0394pj)" ){
  upper_tri <- get_upper_tri(covmat)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)

  plot_cor <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white")

  if (include_values){
    plot_cor <- plot_cor + geom_text(aes(label = round(value, 4)), size = 5, color = "black")
  }

plot_cor <- plot_cor + scale_fill_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      midpoint = 0,
      space = "Lab",
      name = legend_name
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      text = element_text(size = 14, family = "Helvetica")
    ) +
    coord_fixed() +
    labs(x = "", y = "")

  ggsave(output_name,
         plot_cor,
         bg = "white",
         device = cairo_pdf)
}

plot_manhattan_by_site <- function(data, site, today, height, width, quick = T, prefix){

  print(site)

  print("Filter data to one site...")
  data <- data |>
    subset(SITE == site)

  if (quick){
    print("Subset insignificant data...")
    data <- data |>
      subset(pvalue < 0.001)
  }

  print("Compute variables...")
  don <- data %>%
    group_by(CHROM) %>%
    summarise(chr_len=max(POS)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    left_join(data, ., by=c("CHROM"="CHROM")) %>%
    arrange(CHROM, POS) %>%
    mutate(BPcum=POS+tot)

  axisdf = don %>%
    group_by(CHROM) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  print("Plot...")
  ggplot(don, aes(x=BPcum, y=-log10(pvalue))) +
    geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3, shape = 16) +
    geom_point(data=subset(don, highlight=="yes"), color="orange", size=1.5) +
    scale_color_manual(values = rep(c("grey", "black"), 22 )) +
    scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = "-log10(p)", title = paste("Site #", site, sep = " "))

  print("Save png...")
  ggsave(paste("../results/", today, "/", prefix, "_", site, ".jpg", sep = ""), height = height, width = width)
}
