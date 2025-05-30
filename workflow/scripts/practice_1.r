library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)


### <------>### Building Seurat object ###<------>###

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
genome(annotation) <- "hg38"


# load the objects
counts_raw <- Read10X_h5(filename = "data/human_brain_3k_raw_feature_bc_matrix.h5")
# filename = "data/v2/Human Brain Raw Feature Matrix.h5")#human_brain_3k_raw_feature_bc_matrix.h5")
counts_filtered <- Read10X_h5(filename = "data/human_brain_3k_filtered_feature_bc_matrix.h5")
# "../data/v2/human_brain_3k_filtered_feature_bc_matrix.h5")

# visualize some counts_raw rows - sparse matrix format
metadata <- read.csv(
    file = "data/human_brain_3k_per_barcode_metrics.csv",
    header = TRUE,
    row.names = 1
)

# Import the list of cell barcodes that will be used in this notebook
CellSubset <- readRDS("data/lessonCells_byCluster.RDS")

# create a Seurat object containing the RNA adata
multiObj <- CreateSeuratObject(
    counts = counts_raw$`Gene Expression`[, CellSubset],
    assay = "RNA",
    meta.data = metadata,
    project = "Brain",
    min.features = 1,
    min.cells = 1
)

# create ATAC assay and add it to the object
multiObj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts_raw$Peaks[, Cells(multiObj[["RNA"]])], # subset the matrix to keep only cells that are in the RNA object
    sep = c(":", "-"),
    fragments = "data/human_brain_3k_atac_fragments.tsv.gz",
    annotation = annotation,
    genome = "hg38",
    validate.fragments = T,
    min.features = 1,
    min.cells = 1
)

multiObj <- subset(multiObj, cells = Cells(multiObj[["ATAC"]]))
multiObj
