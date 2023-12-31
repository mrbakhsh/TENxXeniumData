% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TENxXeniumData.R
\name{TENxXeniumData}
\alias{TENxXeniumData}
\alias{spe_mouse_brain}
\alias{sfe_mouse_brain}
\alias{spe_human_pancreas}
\alias{sfe_human_pancreas}
\title{TENxXeniumData}
\value{
a \code{\linkS4class{SpatialExperiment}} or
a \code{\linkS4class{SpatialFeatureExperiment}} data objects.
}
\description{
Collection of Xenium spatial transcriptomics datasets from 10X Genomics,
formatted into the Bioconductor classes, the SpatialExperiment or
SpatialFeatureExperiment (SFE). Such Datasets can be used as
examples in packages, tutorials, or for testing purposes.
}
\details{
The following Xenium Spatial Transcriptomics
datasets by 10X Genomics are currently available:
\itemize{
\item{spe_mouse_brain}
\item{sfe_mouse_brain}
\item{spe_human_pancreas}
\item{sfe_human_pancreas}
}
}
\examples{
# initialize hub instance
eh <- ExperimentHub()

# query for TENxXenium datasets
(q <- query(eh, "TENxXenium"))

}
\author{
Matineh Rahmatbakhsh
}
