###############################################################
# Script for creating SpatialExperiments and SpatialFeatureExperiments
# from Xenium data. The constructed data object includes gene expression
# profiles, per-transcript location data, centroid, segmentation
# boundaries (e.g., cell or nucleus boundaries), and image.
###############################################################


######
# Load libraries
######

library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(vroom)
library(dplyr)
library(BumpyMatrix)
library(arrow)
library(sf)
library(tidyr)
library(stringr)
library(RBioFormats)
library(EBImage)
library(DropletUtils)


######
# QCinfo
######

.QCinfo <- function(sfeObject, organism = "human") {
  m <- counts(sfeObject)
  colData(sfeObject)$nCounts <- colSums(m)
  colData(sfeObject)$nGenes <- colSums(m > 0)
  mt_regex <-
    if (organism == "human") "^MT-" else "^Mt-"
  if (any(str_detect(rownames(m), mt_regex))) {
    mito_genes <- str_detect(rownames(m), "^MT-")
    colData(sfeObject)$prop_mito <-
      colSums(m[mito_genes,])/colData(sfeObject)$nCounts
  }
  rowData(sfeObject)$means <- rowMeans(m)
  rowData(sfeObject)$vars <- rowVars(m)
  rowData(sfeObject)$cv2 <-
    rowData(sfeObject)$vars/rowData(sfeObject)$means^2
  sfeObject

}
mols.qv.threshold = 20

###############################################################################

### Mouse Brain: 10x Genomics Xenium In Situ ----------------------------------
# The tiny subset of the raw data was downloaded from 10x website
# https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard
# curl -O https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
# unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip


#### SpatialExperiment (spe) object -------------------------------------------

## open transcripts (i.e., molecules) and select cells that are in cells.csv.gz
cell_info <- vroom(paste0("cells.csv.gz"))

molecule <-
  vroom("transcripts.csv.gz") %>%
  filter(qv > mols.qv.threshold) %>%
  dplyr::select("feature_name","cell_id","x_location","y_location")

diff_colnames <- setdiff(molecule$cell_id, cell_info$cell_id)

if(length(diff_colnames) == 0) {
  molecule
} else {
  molecule <-
    molecule %>% filter(!cell_id %in% diff_colnames )
}

colnames(molecule) <- # rename columns (optional)
  c("gene", "cell", "x", "y")
#molecule <- molecule[1:500,]

mol <- splitAsBumpyMatrix( # convert the transcripts (i.e., molecule) to bumpy matrix
  molecule[, c("x", "y")],
  row = molecule$gene, col = molecule$cell)

## Expression matrix
#sce <- (can you this approach too for experssion matrix)
 # read10xCounts("cell_feature_matrix.h5")
#colnames(sce) <- colData(sce)$Barcode
#sce <- sce[, colnames(sce) %in% molecule$cell]
#sce <- sce[rownames(sce) %in% molecule$gene,]

y <- with(molecule,
          table(gene, cell))

count <- as.matrix(unclass(y))
sce <- SingleCellExperiment(assays = list(counts = count))
counts(sce) <- as(realize(counts(sce)), "dgCMatrix")

## Centroid/spatial information
cell_info <- vroom(paste0("cells.csv.gz")) %>%
  filter(cell_id %in% colnames(sce))

colData(sce) <- cbind(colData(sce), cell_info[,-1])


## Cellular boundary
cell_poly <-
  read_parquet(paste0("cell_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))

cell_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y"))

cells_sf <- cell_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y")) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")


## nuclear boundary
nuc_poly <-
  read_parquet("nucleus_boundaries.parquet") %>%
  filter(cell_id %in% colnames(sce))

nuc_sf <-
  nuc_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y")) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

colData(sce)$cellSeg <- cells_sf$geometry
colData(sce)$nucSeg <- nuc_sf$geometry

spe <- SpatialExperiment(
  assays = list(
    counts = counts(sce),
    molecules = mol),
  colData = colData(sce),
  spatialCoordsNames =
    c("x_centroid", "y_centroid"))

organism = 'mouse'
spe <- .QCinfo(spe, organism = organism)
colData(spe)$total_counts <- NULL

#### SpatialFeatureExperiment (sfe) object ------------------------------------

# cellular boundary
cell_poly <-
  read_parquet(paste0("cell_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))
cells_sf <-
  df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(cells_sf))

## nuclear boundary
nuc_poly <-
  read_parquet(paste0("nucleus_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))
nuc_sf <-
  df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(nuc_sf))


# convert to SpatialFeatureExperiment
sfe <- toSpatialFeatureExperiment(spe)
cellSeg(sfe, withDimnames = FALSE) <- cells_sf
nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
colData(sfe)$cellSeg <- NULL
colData(sfe)$nucSeg <- NULL


#### ADD IMAGE DATA to both data objects -------------------------------------

## ADD IMAGE DATA (ome)
omeImage <-  read.image('morphology.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
image_object <- list()
sample_id <- c()
img_id <- c()
for (channel in 1:num_channels) { # Extract image data from the channel
  img <- SpatialImage(as.raster(imageData(omeImage[,,channel])))
  image_object[[channel]]  <-img
  sample_id[[channel]] <- "sample01" #you can replace sample with ur sampleID
  sample_id <- unlist(sample_id)
  img_id[[channel]] <- paste0("ome", channel) #replace with ur image id
  img_id <- unlist(img_id)

}
scaleFactor = 0.7
scaleFactor = rep(scaleFactor, times = num_channels)
imgData_Ome <- S4Vectors::DataFrame(
  sample_id = sample_id,
  image_id = img_id,
  data = I(image_object),
  scaleFactor = scaleFactor)


## ADD IMAGE DATA (focus.ome)
omeImage <-  read.image('morphology_focus.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
img <-
  SpatialImage(as.raster(imageData(omeImage)))

imgData_focus.ome <- S4Vectors::DataFrame(
  sample_id = "sample01",
  image_id = "focus.ome",
  data = I(list(img)),
  scaleFactor = scaleFactor)


## ADD IMAGE DATA (mip.ome)
omeImage <-  read.image('morphology_mip.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
img <-
  SpatialImage(as.raster(imageData(omeImage)))

imgData_mip.ome <- S4Vectors::DataFrame(
  sample_id = "sample01",
  image_id = "mip.ome",
  data = I(list(img)),
  scaleFactor = scaleFactor)


img_combined <- rbind(imgData_Ome,imgData_focus.ome,imgData_mip.ome)

imgData(spe) <- img_combined # add image to spe
imgData(sfe) <- img_combined # add image to sfe

#### Save datasets ------------------------------------------------------------

saveRDS(spe, file = "spe_mouse_brain.rds")
saveRDS(sfe, file = "sfe_mouse_brain.rds")

###############################################################################


###############################################################################

### Human Pancreas Preview Data: 10x Genomics Xenium In Situ ----------------------------------
# The preview raw data was downloaded from 10x website
# https://www.10xgenomics.com/resources/datasets/human-pancreas-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard
# curl -O https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hPancreas_nondiseased_section/Xenium_V1_hPancreas_nondiseased_section_outs.zip
# unzip Xenium_V1_hPancreas_nondiseased_section_outs.zip


#### SpatialExperiment (spe) object -------------------------------------------
mols.qv.threshold = 20
## open transcripts (i.e., molecules) and select cells that are in cells.csv.gz
cell_info <- vroom(paste0("cells.csv.gz"))

molecule <-
  vroom("transcripts.csv.gz") %>%
  filter(qv > mols.qv.threshold) %>%
  dplyr::select("feature_name","cell_id","x_location","y_location")

diff_colnames <- setdiff(molecule$cell_id, cell_info$cell_id)

if(length(diff_colnames) == 0) {
  molecule
} else {
  molecule <-
    molecule %>% filter(!cell_id %in% diff_colnames )
}

colnames(molecule) <- # rename columns (optional)
  c("gene", "cell", "x", "y")
#molecule <- molecule[1:500,]

mol <- splitAsBumpyMatrix( # convert the transcripts (i.e., molecule) to bumpy matrix
  molecule[, c("x", "y")],
  row = molecule$gene, col = molecule$cell)

## Expression matrix
y <- with(molecule,
          table(gene, cell))
count <- as.matrix(unclass(y))
sce <- SingleCellExperiment(assays = list(counts = count))
counts(sce) <- as(realize(counts(sce)), "dgCMatrix")

## Centroid/spatial information
cell_info <- vroom(paste0("cells.csv.gz")) %>%
  filter(cell_id %in% colnames(sce))

colData(sce) <- cbind(colData(sce), cell_info[,-1])

## Cellular boundary
cell_poly <-
  read_parquet(paste0("cell_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))

cell_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y"))

cells_sf <- cell_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y")) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")


## nuclear boundary
nuc_poly <-
  read_parquet("nucleus_boundaries.parquet") %>%
  filter(cell_id %in% colnames(sce))

nuc_sf <-
  nuc_poly %>%
  st_as_sf(coords = c("vertex_x", "vertex_y")) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

colData(sce)$cellSeg <- cells_sf$geometry
colData(sce)$nucSeg <- nuc_sf$geometry

spe <- SpatialExperiment(
  assays = list(
    counts = counts(sce),
    molecules = mol),
  colData = colData(sce),
  spatialCoordsNames =
    c("x_centroid", "y_centroid"))

organism = 'human'
spe <- .QCinfo(spe, organism = organism)
colData(spe)$total_counts <- NULL

#### SpatialFeatureExperiment (sfe) object ------------------------------------

# cellular boundary
cell_poly <-
  read_parquet(paste0("cell_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))
cells_sf <-
  df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(cells_sf))

## nuclear boundary
nuc_poly <-
  read_parquet(paste0("nucleus_boundaries.parquet")) %>%
  filter(cell_id %in% colnames(sce))
nuc_sf <-
  df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
all(st_is_valid(nuc_sf))


# convert to SpatialFeatureExperiment
sfe <- toSpatialFeatureExperiment(spe)
cellSeg(sfe, withDimnames = FALSE) <- cells_sf
nucSeg(sfe, withDimnames = FALSE) <- nuc_sf
colData(sfe)$cellSeg <- NULL
colData(sfe)$nucSeg <- NULL


#### ADD IMAGE DATA to both data objects -------------------------------------

## ADD IMAGE DATA (ome)
omeImage <-  read.image('morphology.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
image_object <- list()
sample_id <- c()
img_id <- c()
for (channel in 1:num_channels) { # Extract image data from the channel
  img <- SpatialImage(as.raster(imageData(omeImage[,,channel])))
  image_object[[channel]]  <-img
  sample_id[[channel]] <- "sample01" #you can replace sample with ur sampleID
  sample_id <- unlist(sample_id)
  img_id[[channel]] <- paste0("ome", channel) #replace with ur image id
  img_id <- unlist(img_id)
  
}
scaleFactor = 0.7
scaleFactor = rep(scaleFactor, times = num_channels)
imgData_Ome <- S4Vectors::DataFrame(
  sample_id = sample_id,
  image_id = img_id,
  data = I(image_object),
  scaleFactor = scaleFactor)


## ADD IMAGE DATA (focus.ome)
omeImage <-  read.image('morphology_focus.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
img <-
  SpatialImage(as.raster(imageData(omeImage)))

imgData_focus.ome <- S4Vectors::DataFrame(
  sample_id = "sample01",
  image_id = "focus.ome",
  data = I(list(img)),
  scaleFactor = scaleFactor)


## ADD IMAGE DATA (mip.ome)
omeImage <-  read.image('morphology_mip.ome.tif',resolution = 5, series = 1)
num_channels <- omeImage@metadata$coreMetadata$imageCount
num_channels
img <-
  SpatialImage(as.raster(imageData(omeImage)))

imgData_mip.ome <- S4Vectors::DataFrame(
  sample_id = "sample01",
  image_id = "mip.ome",
  data = I(list(img)),
  scaleFactor = scaleFactor)


img_combined <- rbind(imgData_Ome,imgData_focus.ome,imgData_mip.ome)

imgData(spe) <- img_combined # add image to spe
imgData(sfe) <- img_combined # add image to sfe

#### Save datasets ------------------------------------------------------------
saveRDS(spe, file = "spe_human_pancreas.rds")
saveRDS(sfe, file = "sfe_human_pancreas.rds")

###############################################################################
