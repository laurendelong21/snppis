#'This script calculates pathway impact scores of SNPs
#'It takes 3 files as input
#'1. snp.mat : Patients*SNPs (0/1/2) matrix
#'2. pathway_snp_dataframe
#'3. polyphen_sift = polyphen and SIFT scores for each SNPs

library(biomaRt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~Loading all necessary files~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' 1. snp.mat : Patients*SNPs (0/1/2) matrix
load("/home/memon/addneuromed_genetic/processing/4riskmodelvalidation/anm_ctlmci_30k.RData")
snp.mat <- anm_ctrl_30k
rm(anm_ctrl_30k)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#'2. pathway_snp_dataframe
load("/home/memon/addneuromed_genetic/processing/4riskmodelvalidation/pathway_snp_dataframe.RData")
snp2pathway <- pathway_snp_dataframe
snp2pathway <- snp2pathway[,c(2,1)]
names(snp2pathway) <- c("rsID","Pathways")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#'3. polyphen_sift = polyphen and SIFT scores for SNPs in snp.mat using biomart package

if(!file.exists("/home/memon/genetic_analyses/scripts/pathburden/polyphen_sift_anm_30ksnps.RData")){
  cat(sprintf("~~ Retrieving Polyphen and SIFT scores for queried SNPs. ~~\n"))
  rsIDs <- names(snp.mat)[1:100]
  snp.mart = useEnsembl(biomart="snp", dataset="hsapiens_snp")
  polyphen_sift = getBM(attributes=c("refsnp_id", "polyphen_prediction","sift_prediction"), 
                                        mart = snp.mart,filters = "snp_filter", value = rsIDs)
  polyphen_sift$polyphen_prediction <- gsub("^$",NA,polyphen_sift$polyphen_prediction)
  polyphen_sift$sift_prediction <- gsub("^$",NA, polyphen_sift$sift_prediction)
  polyphen_sift = polyphen_sift[which(! is.na(polyphen_sift$sift_prediction)),]
  polyphen_sift = polyphen_sift[which(! is.na(polyphen_sift$polyphen_prediction)),]
  save(polyphen_sift, file="/home/memon/genetic_analyses/scripts/pathburden/polyphen_sift_anm_30ksnps.RData")
} else {
      cat(sprintf("~~'polyphen_sift_anm_30ksnps' already exists, not retrieving again.~~\n"))
      load("/home/memon/genetic_analyses/scripts/pathburden/polyphen_sift_anm_30ksnps.RData")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~Pathway Impact Score calculation function~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pathway.impact.score = function(snp.mat, snp2pathway, polyphen_sift){
  snps.damage = polyphen_sift$refsnp_id[polyphen_sift$polyphen_prediction %in% c("damaging", "possibly damaging", "probably damaging") | polyphen_sift$sift_prediction %in% c("deleterious", "deleterious - low confidence")] # potentially damaging SNPs
  pathways = unique(snp2pathway$Pathway)
  scores = c()
  for(paths in pathways){
    mysnps = intersect(as.character(snp2pathway$rsID[snp2pathway$Pathway == paths]), colnames(snp.mat)) # pathway SNPs in data
    mysnps2 = intersect(mysnps, snps.damage) # pathway SNPs in data that are potentially damaging
    if(length(mysnps2) > 0)
      sc = apply(snp.mat, 1, function(x) sum(x[mysnps2] > 0) / length(mysnps)) # pathway scores for each patient
    else
      sc = rep(0, NROW(snp.mat))
    scores = cbind(scores, sc) 
  }
  colnames(scores) <- pathways
  rownames(scores) <- rownames(snp.mat)
  scores
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~Calcualte Pathway Impact Score~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pathway_burden <- pathway.impact.score(snp.mat, snp2pathway, polyphen_sift)

save(pathway_burden, file = "/home/memon/addneuromed_genetic/processing/4riskmodelvalidation/anm_ctlmci_30K_pathway_score.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
