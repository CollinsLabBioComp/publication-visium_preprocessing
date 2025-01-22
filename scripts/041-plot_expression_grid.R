library(optparse)
library(tidyverse)
library(cowplot)
library(grid)

option_list = list(
  make_option(c("-p", "--path"), action = "store", default = "", type = 'character',
              help = "Path to spaceranger counts outs"),
  make_option(c("-r", "--radii"), action = "store", default = "", type = 'character',
              help = "Decontamination radii"),
  make_option(c("-g", "--gene"), action = "store", default = "", type = 'character',
              help = "Gene expression to plot"),
  make_option(c("-k","--kernel"), action = "store", default = "", type = 'character',
              help = "Kernels to mdoel decontamination")
  )
opt = parse_args(OptionParser(option_list = option_list))

print(opt$p)
print(opt$r)
print(opt$g)

cwd <- getwd()
print(cwd)
radii <- opt$r

split <- strsplit(radii, ' ')

kernel <- opt$k
gene <- opt$g

PATH <- file.path(cwd, opt$p)
path_split <- strsplit(PATH, '/')
sample <- path_split[[1]][9]
print(PATH)

plots <- list()

spaceranger <- read_csv(paste0(PATH, 'spot_', gene, '_expression.csv'))
spaceranger

plots[[1]] <- ggplot(spaceranger, aes(x,y, fill=expr, color=expr)) +
    geom_point(size=1.5) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    labs(x='X', y='Y', color='log(counts+1)',fill='log(counts+1)',title=paste0(gene, ' Expression Space Ranger')) +
    theme_bw()

spaceranger_norm <- read_csv(paste0(PATH, 'spot_', gene, '_norm_expression.csv'))

plots[[2]] <- ggplot(spaceranger_norm, aes(x,y, fill=expr, color=expr)) +
    geom_point(size=1.5) +
    scale_fill_viridis_c() +
    scale_color_viridis_c() +
    labs(x='X', y='Y', color='log(cp10k+1)',fill='log(cp10k+1)',title=paste0('Normalized ', gene, ' Expression Space Ranger')) +
    theme_bw()

for (i in split[[1]]) {
    i <- as.numeric(i)
    spotclean <- read_csv(paste0(PATH, 'radius=', i, '/spot_', gene, '_expression.csv'))
    index = which(split[[1]] == i) + 2
    print(index)
    print(spotclean)

    plots[[index]] <- ggplot(spotclean, aes(x,y, fill=expr, color=expr)) +
        geom_point(size=1.5) +
        scale_fill_viridis_c() +
        scale_color_viridis_c() +
        labs(x='X',
            y='Y',
            color='log(sc_counts+1)',
            fill='log(sc_counts+1)',
            title=paste0('SpotCleaned radius=', i, ' ', gene, ' Expression')) +
        theme_bw() +
        theme(legend.title = element_text(size = 10))
}

title <- ggdraw() + 
  draw_label(
    paste0(sample, " SpotClean v1.4.1 kernel=",kernel," max_iterations=30"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

plot_grid_matrix <- plot_grid(plotlist = plots, nrow = 2)

plot_grid_matrix_title <- plot_grid(title, plot_grid_matrix, ncol = 1, rel_heights = c(0.1,1))


ggsave(plot_grid_matrix_title, file = paste0(PATH, 'combined_', gene, '_expression.png'), width = 22, height = 10)
