library(ggplot2)
library(stringr)
library(grid)

myplot <- function(df, 
                   title = "Pathway Analysis",
                   save_prefix = "kegg_plot",
                   output_dir = ".",
                   width = 9,
                   height = 6) {
  

  df$Description <- str_wrap(df$Description, width = 30)
  

  df$GeneRatio <- as.numeric(str_split_fixed(as.character(df$GeneRatio), "/", 2)[,1]) /
                as.numeric(str_split_fixed(as.character(df$GeneRatio), "/", 2)[,2])
  

  count_breaks<- c(
    min(df$Count, na.rm = TRUE),
    round(mean(df$Count, na.rm = TRUE)),
    max(df$Count, na.rm = TRUE)
)
  
  pvalue_min <- min(df$pvalue, na.rm = TRUE)
  pvalue_max <- max(df$pvalue, na.rm = TRUE)
  pvalue_middle <- (pvalue_min + pvalue_max) / 2
  
  color_breaks <- c(pvalue_min, pvalue_middle, pvalue_max)
  color_values <- scales::rescale(color_breaks) 
  color_labels <- sprintf("%.2e", color_breaks)
  
  p <- ggplot(df, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradientn(
      colours = c("#DC0000FF", "#F79D1E99", "#357EBD99"),
      values = color_values,
      breaks = color_breaks,
      labels = color_labels,
      limits = c(pvalue_min, pvalue_max),
      guide = guide_colorbar(
        reverse = TRUE,
        title = "P value",
        order = 1,
        barheight = unit(4, "cm"),
        barwidth = unit(0.8, "cm"))
    ) +
    scale_size_continuous(
        range = c(3, 8),
        breaks = count_breaks,
        labels = as.character(count_breaks),
        guide = guide_legend(
            title = "Count",
            override.aes = list(color = "black")
        )
    ) +
    labs(title = title, x = "Gene Ratio", y = "Descriptions") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(margin = margin(r = 10))
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  ggsave(
    file.path(output_dir, paste0(save_prefix, ".png")),
    plot = p,
    width = width,
    height = height
  )
  ggsave(
    file.path(output_dir, paste0(save_prefix, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
  return(df)
}
