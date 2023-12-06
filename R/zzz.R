#' @importFrom utils read.csv
#' @importFrom ExperimentHub createHubAccessors
#' @import SpatialExperiment
#' @import SpatialFeatureExperiment

.onLoad <- function(libname, pkgname) {
  fl <- system.file("extdata", "metadata.csv", package=pkgname)
  titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
  createHubAccessors(pkgname, titles)
}
