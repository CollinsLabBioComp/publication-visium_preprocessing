library(optparse)
library(tidyverse)
library(cowplot)
library(grid)

option_list = list(
  make_option(c("-p", "--path"), action = "store", default = "", type = 'character',
              help = "Output directory for plots"),
  make_option(c("-g", "--gene"), action = "store", default = "", type = 'character',
              help = "Gene expression to plot")
  )
opt = parse_args(OptionParser(option_list = option_list))

print(opt$p)
print(opt$g)

cwd <- getwd()
print(cwd)

gene <- opt$g

PATH <- file.path(cwd, opt$p)
path_split <- strsplit(PATH, '/')
sample <- path_split[[1]][9]
print(PATH)

plots <- list()

spaceranger <- read_csv(paste0(PATH, 'spot_', gene, '_expression.csv'))
spaceranger

plots[[1]] <- ggplot(spaceranger, aes(x,y, fill=expr, color=expr)) +
    geom_point(size=2.5) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    labs(x='X', y='Y', color='log(counts+1)',fill='log(counts+1)',title=paste0(gene, ' Expression Space Ranger')) +
    theme_bw()

spaceranger_norm <- read_csv(paste0(PATH, 'spot_', gene, '_norm_expression.csv'))

plots[[2]] <- ggplot(spaceranger_norm, aes(x,y, fill=expr, color=expr)) +
    geom_point(size=2.5) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    labs(x='X', y='Y', color='log(cp10k+1)',fill='log(cp10k+1)',title=paste0('Normalized ', gene, ' Expression Space Ranger')) +
    theme_bw()


title <- ggdraw() + 
  draw_label(
    paste0(sample, " Normalization"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

plot_grid_matrix <- plot_grid(plotlist = plots, nrow = 1)

plot_grid_matrix_title <- plot_grid(title, plot_grid_matrix, ncol = 1, rel_heights = c(0.1,1))

ggsave(plots[[1]], file = paste0(PATH, gene, '_expression.png'))
ggsave(plots[[2]], file = paste0(PATH, 'normalized_', gene, '_expression.png'))
ggsave(plot_grid_matrix_title, file = paste0(PATH, 'compare_normalization_', gene, '.png'), width = 16, height = 8)
