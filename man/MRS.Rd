% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRS.R
\name{MRS}
\alias{MRS}
\title{Microbial risk score framework}
\usage{
MRS(discovery.data, validation.data, GroupID, DA.method, measurement)
}
\arguments{
\item{discovery.data}{A phyloseq-class object, which consists of a feature table (observed count table), a sample metadata,
a taxonomy table (optional), and a phylogenetic tree (optional). See \code{\link[phyloseq]{phyloseq}} for more details.}

\item{validation.data}{A phyloseq-class object, which consists of a feature table (observed count table or relative abundance table),
a sample metadata, a taxonomy table (optional), and a phylogenetic tree (optional).
See \code{\link[phyloseq]{phyloseq}} for more details.}

\item{GroupID}{The name of the group variable of interest in metadata.}

\item{DA.method}{The microbial differential abundance method which is employed to identify the candidate taxa for MRS
construction based on the discovery dataset (the first step of the MRS framework).
Three top-performing methods \href{https://www.nature.com/articles/s41467-020-17041-7}{ANCOMBC},
\href{https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html}{ALDEx2},
and \href{https://huttenhower.sph.harvard.edu/maaslin/}{Maaslin2} are available: c("ancombc", "ALDEx2","Maaslin2").}

\item{measurement}{The diverity index which is employed for MRS calculation.
Three widely used indices are available: c("shannon","simpson","observed").}
}
\value{
A list with components:
\itemize{
\item{ Cutoff: An optimal p-value cutoff (a numeric value) for identifying the taxa used for construction of MRS.}
\item{ Taxa used for MRS: The specific taxa used for construction of MRS.}
\item{ AUC: A numeric matrix which reports the AUC values and 95\% confidence intervals in
the discovery and validation cohorts, respectively.}}
}
\description{
\code{MRS} Construction and evaluation of MRS.
}
\details{
This function identifys the candidate taxa used for constructing MRS and calculates MRS in the discovery cohort,
and evaluates MRS in both discovery and validation cohorts in terms of AUC.
If there is no validation cohort, cross-validation can be performed at the pre-processing step to have both cohorts.

\strong{Reference}: Wang C, Segal L, Hu J, Zhou B, Hayes R, Ahn J, and Li H (2022).
Microbial risk score: capturing microbial characteristics, integrating multi-omics profiling,
and predicting disease risk. \emph{Microbiome}.
}
\examples{

# require(phyloseq)

## Evaluation of MRS in terms of comparison between Healthy and Nonhealthy ##
## using ANCOMBC method and Shannon index
#  discovery=GMHI[[1]];
# validation=GMHI[[2]];

# res=MRS(discovery, validation, GroupID="Group", DA.method="ancombc", measurement="shannon")
#  AUC=res[[3]]

 ## using ALDEx2 method and Shannon index

# res=MRS(discovery, validation, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
# AUC=res[[3]]

## Evaluation of MRS in terms of comparison between Healthy and a specific disease ##
## Healthy vs. CA

# discovery.sub=prune_samples(sample_data(discovery)$Group1 \%in\% c("Healthy","CA"),discovery)
# validation.sub=prune_samples(sample_data(validation)$Group1 \%in\% c("Healthy","CA"),validation)

# res=MRS(discovery.sub, validation.sub, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
# AUC=res[[3]]

}
