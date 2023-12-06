# ------------------------------
# Generate metadata spreadsheet
# -----------------------------

# metadata for all datasets

df_all <- data.frame(
  BiocVersion = "3.18",
  Genome = NA,
  SourceVersion = NA,
  Coordinate_1_based = NA,
  DataProvider = NA,
  Maintainer = "Matineh Rahmatbakhsh <matinerb.94@gmail.com>",
  DispatchClass = "Rds",
  stringsAsFactors = FALSE
)

df_spe_mouse_brain <- cbind(
  df_all,
  Title = "spe_mouse_brain",
  Description = "Example Xenium Tiny subset dataset from 10x Genomics provided in the Fresh Frozen Mouse Brain for Xenium Explorer Demo",
  SourceUrl = "https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  RDataPath = "TENxXeniumData/spe_mouse_brain.rds",
  RDataClass = "SpatialExperiment",
  SourceType = "Zip",
  stringsAsFactors = FALSE
)


df_sfe_mouse_brain <- cbind(
  df_all,
  Title = "sfe_mouse_brain",
  Description = "Example Xenium Tiny subset dataset from 10x Genomics provided in the Fresh Frozen Mouse Brain for Xenium Explorer Demo",
  SourceUrl = "https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard",
  Species = "Mus musculus",
  TaxonomyId = "10090",
  RDataPath = "TENxXeniumData/sfe_mouse_brain.rds",
  RDataClass = "SpatialFeatureExperiment",
  SourceType = "Zip",
  stringsAsFactors = FALSE
)



# combine and save as .csv spreadsheet file

df_combined <- rbind(
  df_spe_mouse_brain,
  df_sfe_mouse_brain
)

write.csv(df_combined, file = "metadata.csv", row.names = FALSE)
