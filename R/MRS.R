#' Microbial risk score framework
#'
#' \code{MRS} Construction and evaluation of MRS.
#'
#' This function identifys the candidate taxa used for constructing MRS and calculates MRS in the discovery cohort,
#' and evaluates MRS in both discovery and validation cohorts in terms of AUC.
#' If there is no validation cohort, cross-validation can be performed at the pre-processing step to have both cohorts.
#'
#'\strong{Reference}: Wang C, Segal L, Hu J, Zhou B, Hayes R, Ahn J, and Li H (2022).
#'Microbial risk score: capturing microbial characteristics, integrating multi-omics profiling,
#'and predicting disease risk. \emph{Microbiome}.
#'
#' @param discovery.data A phyloseq-class object, which consists of a feature table (observed count table), a sample metadata,
#' a taxonomy table (optional), and a phylogenetic tree (optional). See \code{\link[phyloseq]{phyloseq}} for more details.
#'
#' @param validation.data A phyloseq-class object, which consists of a feature table (observed count table or relative abundance table),
#' a sample metadata, a taxonomy table (optional), and a phylogenetic tree (optional).
#' See \code{\link[phyloseq]{phyloseq}} for more details.
#'
#' @param GroupID The name of the group variable of interest in metadata.
#'
#' @param DA.method The microbial differential abundance method which is employed to identify the candidate taxa for MRS
#' construction based on the discovery dataset (the first step of the MRS framework).
#' Three top-performing methods \href{https://www.nature.com/articles/s41467-020-17041-7}{ANCOMBC},
#'\href{https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html}{ALDEx2},
#'and \href{https://huttenhower.sph.harvard.edu/maaslin/}{Maaslin2} are available: c("ancombc", "ALDEx2","Maaslin2").
#'
#' @param measurement The diverity index which is employed for MRS calculation.
#' Three widely used indices are available: c("shannon","simpson","observed").
#'
#' @return A list with components:
#' \itemize{
#'  \item{ Cutoff: An optimal p-value cutoff (a numeric value) for identifying the taxa used for construction of MRS.}
#' \item{ Taxa used for MRS: The specific taxa used for construction of MRS.}
#'\item{ AUC: A numeric matrix which reports the AUC values and 95% confidence intervals in
#'the discovery and validation cohorts, respectively.}}
#'
#' @examples
#'
#' require(phyloseq)
#'
#' ## Evaluation of MRS in terms of comparison between Healthy and Nonhealthy ##
#' ## using ANCOMBC method and Shannon index
#'  discovery=GMHI[[1]];
#'
#'  validation=GMHI[[2]];
#'
#'  res=MRS(discovery, validation, GroupID="Group", DA.method="ancombc", measurement="shannon")
#'
#'  AUC=res[[3]]
#'
#'  ## using ALDEx2 method and Shannon index
#'
#'  res=MRS(discovery, validation, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
#'
#'  AUC=res[[3]]
#'
#' ## Evaluation of MRS in terms of comparison between Healthy and a specific disease ##
#' ## Healthy vs. CA
#'
#'  discovery.sub=prune_samples(sample_data(discovery)$Group1 %in% c("Healthy","CA"),discovery)
#'  validation.sub=prune_samples(sample_data(validation)$Group1 %in% c("Healthy","CA"),validation)
#'
#'  res=MRS(discovery.sub, validation.sub, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
#'  AUC=res[[3]]
#'
#'@importFrom magrittr "%>%"
#'@import ANCOMBC
#'@import Maaslin2
#'@import ALDEx2
#'@import phyloseq
#'@import pROC
#'@import vegan
#'
#' @export
MRS=function(discovery.data, validation.data, GroupID,
             DA.method,measurement) {

  if(DA.method=="ancombc") {

  out = ancombc2(data = discovery.data, fix_formula = GroupID,
                p_adj_method = "BH",  #zero_cut = 0.90, #lib_cut = 1000,
                group = GroupID, struc_zero = TRUE, neg_lb = TRUE,
                alpha = 0.05, global = FALSE)

  res0 = out$res
  res = res0[,grep("p_",colnames(res0))[2]] %>% as.matrix() %>% data.frame()

  colnames(res)="pvalue"
  res$prank=rank(res$pvalue,ties.method ="min")
  res$taxon=rownames(res0)}

  if(DA.method=="ALDEx2") {

    OTU=otu_table(discovery.data)
    meta=sample_data(discovery.data)

    RA_clr=aldex.clr(reads=t(OTU),conds = as.character(unlist(meta[,GroupID])))
    res=aldex.ttest(RA_clr)


    res = res %>% as.matrix() %>% data.frame()
    colnames(res)[1]="pvalue"
    res$prank=rank(res$pvalue,ties.method ="min")
    res$taxon=rownames(res) }


  if(DA.method=="Maaslin2") {

    OTU=data.frame(otu_table(discovery.data))
    meta=data.frame(sample_data(discovery.data))
    rownames(OTU)=rownames(meta)

    res=Maaslin2(
      input_data = OTU,
      input_metadata = meta,
      output = "demo_output", transform="AST",
      fixed_effects = GroupID,
      plot_heatmap =FALSE,
      plot_scatter = FALSE)


    res=res$results[,c("feature","pval")]
    res = res %>% as.matrix() %>% data.frame()
    colnames(res)[2]="pvalue"
    res$prank=rank(res$pvalue,ties.method ="min")
    colnames(res)[1]="taxon" }


  AUC.res=NULL

  for(nTop in  seq(5, nrow(res),1)) {

    community.used=prune_taxa(res$taxon[res$prank<=nTop], discovery.data)

    AA=otu_table(community.used)

    if(measurement!="observed") {

      erDF = diversity(AA,index=measurement)
      Value=as.numeric(erDF)} else Value=as.numeric(apply(otu_table(community.used),1, function(x) sum(x>0)))

    roc1<- roc(status~Value,
               data=data.frame(status=as.character(unlist(sample_data(discovery.data)[,GroupID])),
                                            Value=Value))
    AUC.res=c(AUC.res,roc1$auc)
  }


    pvalue.opt=seq(5, nrow(res),1)[which.max(AUC.res)]

    Taxa_identified=res$taxon[res$prank<=pvalue.opt]
    community.used=prune_taxa(Taxa_identified, discovery.data)
    AA=otu_table(community.used)

    if(measurement!="observed") {

      erDF = diversity(AA,index=measurement)
      Value=as.numeric(erDF) } else Value=as.numeric(apply(otu_table(community.used),1, function(x) sum(x>0)))

    roc1<- roc(status~Value,
               data=data.frame(status=as.character(unlist(sample_data(discovery.data)[,GroupID])),
                               Value=Value))
    aa0=ci(roc1)


    community.used=prune_taxa(Taxa_identified, validation.data)
    AA=otu_table(community.used)

    if(measurement!="observed") {

      erDF = diversity(AA,index=measurement)
      Value=as.numeric(erDF) } else Value=as.numeric(apply(otu_table(community.used),1, function(x) sum(x>0)))

    roc1<- roc(status~Value,
               data=data.frame(status=as.character(unlist(sample_data(validation.data)[,GroupID])),
                               Value=Value))

    aa1=ci(roc1)

    AUC=round(rbind(c(aa0[2],aa0[1],aa0[3]),
          c(aa1[2],aa1[1],aa1[3])),3)

    rownames(AUC)=c("Discovery cohort","Validation cohort")
    colnames(AUC)=c("AUC","LCI","UCI")


    pcutoff=res$pvalue[res$prank==pvalue.opt]

    Allresults=list(pcutoff,Taxa_identified,AUC)
    names(Allresults)=c("Cutoff","Taxa used for MRS","AUC")

    return(Allresults)

  }
