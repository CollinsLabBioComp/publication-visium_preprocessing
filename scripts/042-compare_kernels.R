library(optparse)
library(tidyverse)
library(cowplot)
library(grid)

option_list = list(
  make_option(c("-p", "--path"), action = "store", default = "", type = 'character',
              help = "Path to spaceranger counts outs"),
  make_option(c("-g", "--gene"), action = "store", default = "", type = 'character',
              help = "Gene expression to plot"),
  make_option(c("-k","--kernels"), action = "store", default = "", type = 'character',
              help = "Kernels to model decontamination")
  )
opt = parse_args(OptionParser(option_list = option_list))

print(opt$p)
print(opt$g)

cwd <- getwd()
print(cwd)

kernels <- opt$k
kernels_split <- strsplit(kernels, ' ')
print(kernels)


gene <- opt$g

PATH <- file.path(cwd, opt$p)
path_split <- strsplit(PATH, '/')
sample <- path_split[[1]][9]
print(PATH)

plots <- list()

for (i in kernels_split[[1]]) {
    spotclean <- read_csv(paste0(PATH, 'kernel=', i, '/radius=10/spot_', gene, '_expression.csv'))
    index = which(kernels_split[[1]] == i)
    print(index)
    print(spotclean)

    plots[[index]] <- ggplot(spotclean, aes(x,y, fill=expr, color=expr)) +
        geom_point(size=1.5) +
        scale_fill_viridis_c() +
        scale_color_viridis_c() +
        labs(x='X',
            y='Y',
            color='log(cp10k+1)',
            fill='log(cp10k+1)',
            title=paste0('SpotCleaned radius=10 kernel=', i, ' ', gene)) +
        theme_bw()
}

for (i in kernels_split[[1]]) {
    spotclean <- read_csv(paste0(PATH, 'kernel=', i, '/radius=30/spot_', gene, '_expression.csv'))
    index = which(kernels_split[[1]] == i) + 4
    print(index)
    print(spotclean)

    plots[[index]] <- ggplot(spotclean, aes(x,y, fill=expr, color=expr)) +
        geom_point(size=1.5) +
        scale_fill_viridis_c() +
        scale_color_viridis_c() +
        labs(x='X',
            y='Y',
            color='log(cp10k+1)',
            fill='log(cp10k+1)',
            title=paste0('SpotCleaned radius=30 kernel=', i, ' ', gene)) +
        theme_bw()
}

title <- ggdraw() + 
  draw_label(
    paste0(sample, " SpotClean v1.4.1 Kernel Comparison max_iterations=30"),
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

ggsave(plot_grid_matrix_title, file = paste0(PATH, 'kernels_', gene, '_expression.png'), width = 20, height = 10)
