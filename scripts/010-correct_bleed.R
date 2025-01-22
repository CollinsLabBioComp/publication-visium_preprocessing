library(Seurat)
library(SpotClean)
library(SeuratDisk)
library(optparse)
library(sceasy)
library(cowplot)
library(S4Vectors)
library(ggplot2)
library(grid)

option_list = list(
  make_option(c("-i", "--input"), action = "store", default = "", type = 'character',
              help = "Path to spaceranger counts outs"),
  make_option(c("-o", "--output"), action = "store", default = "", type = 'character',
              help = "Output directory"),
  make_option(c("-r", "--radius"), action = "store", default = "", type = 'integer',
              help = "Decontamination radius"),
  make_option(c("-k", "--kernel"), action = "store", default = "", type = 'character',
              help = "Kernel to model contamination")
  )
opt = parse_args(OptionParser(option_list = option_list))

# =====================================
# Run
# =====================================

cwd <- getwd()

if (grepl("HD",opt$i)) {
  IN_PATH = file.path(cwd, opt$i, "/binned_outputs/square_002um")
  
  tissue_parquet_file <- read_parquet(file = paste0(IN_PATH, "/spatial/tissue_positions.parquet"))
  write.csv(tissue_parquet_file, 
            file = paste0(IN_PATH, "/spatial/tissue_positions.csv"), 
            row.names = FALSE,
            quote = FALSE)

  slide_info <- read10xSlide(tissue_csv_file = paste0(IN_PATH, "/spatial/tissue_positions.csv"),
                             tissue_img_file = paste0(IN_PATH, "/spatial/tissue_lowres_image.png"),
                             scale_factor_file = paste0(IN_PATH, "/spatial/scalefactors_json.json"))
} else {
  IN_PATH = file.path(cwd, opt$i)
  slide_info <- read10xSlide(tissue_csv_file = paste0(IN_PATH, "/spatial/tissue_positions.csv"),
                             tissue_img_file = paste0(IN_PATH, "/spatial/tissue_lowres_image.png"),
                             scale_factor_file = paste0(IN_PATH, "/spatial/scalefactors_json.json"))
}

radius <- opt$r
OUT_PATH = file.path(cwd, opt$o)
split <- strsplit(IN_PATH, '/')
sample <- split[[1]][9]
kernel <- opt$k

# Read in raw count matrix and slide information
raw_mat <- read10xRaw(paste0(IN_PATH,'/raw_feature_bc_matrix'))

slide_obj <- createSlide(count_mat = raw_mat, 
                        slide_info = slide_info,
                        gene_cutoff = 0.01)


# Create SpotClean intermediate plots
plots <- list()
plots[[1]] <- visualizeSlide(slide_obj = slide_obj, title=paste0(sample, ' Slide Image'))
plots[[2]] <- visualizeLabel(slide_obj,"tissue", title=paste0(sample, ' In-Tissue Barcodes'))
metadata(slide_obj)$slide$total_counts <- Matrix::colSums(raw_mat)
plots[[3]] <- visualizeHeatmap(slide_obj, "total_counts", title=paste0(sample, ' Total Counts'))
plots[[4]] <- visualizeHeatmap(slide_obj, 'INS', title=paste0(sample, ' INS Expression'))

decont_obj <- spotclean(slide_obj, candidate_radius = radius, kernel = kernel)

plots[[5]] <- visualizeHeatmap(decont_obj, 'INS', title=paste0(sample, ' Cleaned INS Expression'))


# Plot grid and save
plot_grid_matrix <- plot_grid(plotlist = plots, nrow = 1)
ggsave(plot_grid_matrix, file = paste0(OUT_PATH, '/kernel=', kernel, '/', '/radius=', radius, '/', sample, '_SpotClean_timeline.png'), width = 40, height = 8)

# =====================================
# Define functions
# =====================================

export_decont <- function(decont_obj, spatial_folder, out_folder, file_name) {
  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  seurat_obj <- convertToSeurat(decont_obj, spatial_folder)
  seurat_obj[['Spatial']] <- as(object = seurat_obj[['Spatial']], Class = "Assay")

  updated_seurat_obj <- UpdateSeuratObject(seurat_obj)
  sceasy::convertFormat(updated_seurat_obj, from="seurat", to="anndata",
                       outFile=paste0(out_folder, file_name, ".h5ad"), assay = 'Spatial')
}

export_decont(decont_obj, paste0(IN_PATH, "/spatial/"), paste0(OUT_PATH, "/kernel=", kernel, "/radius=", radius, "/"),  "cleaned_feature_bc_matrix")
