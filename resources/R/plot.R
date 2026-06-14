box::use(ggplot2[...])
box::use(reshape2[melt])
box::use(dplyr[...])
box::use(grDevices[...])
box::use(scico[...])
box::use(ggpubr[...])

#' Title
#'
#' @param mat 
#'
#' @return
#' @export
#'
#' @examples
get_upper_tri <- function(mat){
  mat[lower.tri(mat)]<- NA
  diag(mat) <- NA
  return(mat)
}


#' Title
#'
#' @param melted_cormat 
#' @param Var2 
#' @param Var1 
#' @param value 
#' @param include_values 
#' @param legend_name 
#' @param lwr 
#' @param upr 
#' @param scico_palette 
#' @param midpoint 
#' @param direction 
#' @param height 
#' @param width 
#' @param output_name 
#'
#' @return
#' @export
#'
#' @examples
plot_pairs_heatmap <- function(melted_cormat,
                               Var2,
                               Var1,
                               value,
                               include_values = F,
                               lwr,
                               upr,
                               scico_palette,
                               midpoint = 0,
                               direction = 1,
                               # low_col = "red",
                               # mid_col = "black",
                               # high_col = "blue",
                               height,
                               width,
                               legend_name,
                               output_name) {
  
  plot_cor <- ggplot(data = melted_cormat, 
                     aes({{ Var2 }}, {{ Var1 }}, fill = {{ value }})) +
    geom_tile(color = "white")
  
  plot_cor <- plot_cor + scale_fill_scico(
    # low = low_col,
    # mid = mid_col,
    # high = high_col,
    palette = scico_palette,
    midpoint = midpoint,
    direction = direction,
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
  
  if(include_values){
    plot_cor <- plot_cor + geom_text(aes(label = if_else((lwr < 0) & (upr > 0), "", "*")), 
                                     color = "white", 
                                     size = 8, 
                                     vjust = 0.8) 
  }
  
  ggsave(output_name,
         plot_cor,
         height = height,
         width = width,
         bg = "white")
}

#' Plot heatmap from var-cov matrix
#'
#' @param covmat 
#' @param output_name 
#' @param include_values 
#' @param low_col 
#' @param mid_col 
#' @param high_col 
#' @param legend_name 
#'
#' @return
#' @export
#'
#' @examples
plot_var_cov_matrix <- function(covmat,
                                output_name,
                                include_values = TRUE,
                                low_col = "#65014B",
                                mid_col = "#F5F0F0",
                                high_col = "#0C4C00",
                                legend_name = "Cov(\u0394pi, \u0394pj)",
                                include_tick_labels = TRUE,
                                width = 7,
                                height = 7,
                                tile_or_raster = "tile",
                                discrete_or_continuous = "continuous") {
  
  upper_tri <- get_upper_tri(covmat)
  
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  plot_cor <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))
  
  if (tile_or_raster == "tile") {
    plot_cor <- plot_cor + geom_tile(color = "white")
  } else if (tile_or_raster == "raster") {
    plot_cor <- plot_cor + geom_raster()
  } else {
    stop("tile_or_raster invalid.")
  }
  
  if (include_values) {
    plot_cor <- plot_cor + geom_text(aes(label = round(value, 4)),
                                     size = 5,
                                     color = "black")
  }
  
  if (discrete_or_continuous == "discrete") {
    plot_cor <- plot_cor + scale_fill_scico_d(
      palette = "bam",
      begin = 0.1,
      end = 0.9,
      name = legend_name
    )
  } else if (discrete_or_continuous == "continuous") {
    plot_cor <- plot_cor + scale_fill_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      midpoint = 0,
      space = "Lab",
      name = legend_name
    )
  } else {
    stop("discrete_or_continuous invalid.")
  }
  
  plot_cor <- plot_cor + theme_minimal()
  
  if (include_tick_labels) {
    plot_cor <- plot_cor + theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      text = element_text(size = 14, family = "Helvetica")
    )
  } else {
    plot_cor <- plot_cor +
      theme(axis.text.x = element_blank(), 
            axis.text.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank())
  }
  
  # remove axis names
  plot_cor <- plot_cor + coord_fixed() +
    labs(x = "", y = "")
  
  ggsave(
    output_name,
    plot_cor,
    bg = "white",
    device = cairo_pdf,
    width = width,
    height = height
  )
}

#' Manhattan plot from snp table of stats
#'
#' @param data 
#' @param site 
#' @param today 
#' @param height 
#' @param width 
#' @param quick 
#' @param prefix 
#'
#' @return
#' @export
#'
#' @examples
plot_manhattan_by_site <- function(data, site, today, height, width, quick = F, save_pdf = T, prefix){

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
    geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3, shape = 16) +
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
  ggsave(paste("../results/by-date/", today, "/", prefix, "_", site, ".jpg", sep = ""), height = height, width = width)

  if(save_pdf){
    print("Save pdf...")
    ggsave(paste("../results/by-date/", today, "/", prefix, "_", site, ".pdf", sep = ""), height = height, width = width)
  }
}
