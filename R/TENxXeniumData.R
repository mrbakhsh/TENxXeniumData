  #' @name TENxXeniumData
  #' @title TENxXeniumData
  #' @aliases
  #' spe_mouse_brain
  #' sfe_mouse_brain
  #' @description
  #' Collection of Xenium spatial transcriptomics datasets from 10X Genomics,
  #' formatted into the Bioconductor classes, the SpatialExperiment or
  #' SpatialFeatureExperiment (SFE). Such Datasets can be used as
  #' examples in packages, tutorials, or for testing purposes.
  #' @details
  #' The following Xenium Spatial Transcriptomics
  #' datasets by 10X Genomics are currently available:
  #' \itemize{
  #' \item{spe_mouse_brain}
  #' \item{sfe_mouse_brain}
  #' }
  #' @return
  #' a \code{\linkS4class{SpatialExperiment}} or
  #' a \code{\linkS4class{SpatialFeatureExperiment}} data objects.
  #'
  #' @examples
  #' # initialize hub instance
  #' eh <- ExperimentHub()
  #'
  #' # query for TENxXenium datasets
  #' (q <- query(eh, "TENxXenium"))
  #'
  #' @author Matineh Rahmatbakhsh
  NULL



