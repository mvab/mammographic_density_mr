
kable_it<-function(df){
  library(kableExtra)
  df %>% 
    tidy_pvals %>% 
    kable(.) %>%
    kable_styling()
}
#dat %>% kable_it()


tidy_pvals<-function(df){
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}


mr_scatter_plot_manual <- function (mr_results, dat) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome -  se.outcome, 
                                                               ymax = beta.outcome + se.outcome), 
                                                               colour = "grey", width = 0) +
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure -  se.exposure,
                                                                xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
                           ggplot2::geom_point() + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                           slope = b, colour = method), show.legend = TRUE) + 
                           ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                   "#6a3d9a", "#ffff99", "#b15928")) +
                           ggplot2::labs(colour = "MR Test", 
                                          x = paste("SNP effect on", d$exposure[1]),
                                          y = paste("SNP effect on",  d$outcome[1])) + 
                           
                           ggplot2::theme(legend.position = "top",   legend.direction = "vertical") + 
                          
                           
                           theme_bw() + 
                           ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
                       })
  mrres
}


mr_leaveoneout_plot_manual <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) +
                                                ggplot2::geom_vline(xintercept = 0,  linetype = "dotted") +
                                                ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo,  xmax = up, 
                                                                                     size = as.factor(tot), 
                                                                                     colour = as.factor(tot)),
                                                                        height = 0) +
                                                ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in%  "")), colour = "grey") + 
                                                
                                                ggplot2::scale_colour_manual(values = c("black",  "red")) + 
                                                ggplot2::scale_size_manual(values = c(0.3,   1)) +
                                                theme_bw()+
                                                ggplot2::theme(legend.position = "none", 
                                                               axis.text.y = ggplot2::element_text(size = 8), 
                                                               axis.ticks.y = ggplot2::element_line(size = 0), 
                                                               axis.title.x = ggplot2::element_text(size = 8)) + 
                                                ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", 
                                                                                 d$exposure[1], "' on '", d$outcome[1], "'"))
                                            })
  res
}



mr_forest_plot_outliers <- function (d,  outliers_list, outliers_colour, method) {
    # this is a modified version of 2SMR function 

     requireNamespace("ggplot2", quietly = TRUE)
     requireNamespace("plyr", quietly = TRUE)

     if (sum(!grepl("All|Outlier", d$SNP)) < 2) {
       return(blank_plot("Insufficient number of SNPs"))}

     levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
     levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
     am <- grep("All|Outlier", d$SNP, value = TRUE)
     d$up <- d$b + 1.96 * d$se
     d$lo <- d$b - 1.96 * d$se
     d$tot <- 0.01
     d$tot[d$SNP %in% am] <- 1
     d$SNP <- as.character(d$SNP)
     nom <- d$SNP[!d$SNP %in% am]
     nom <- nom[order(d$b)]
     d <- rbind(d, d[nrow(d), ])
     d$SNP[nrow(d) - 1] <- ""
     d$b[nrow(d) - 1] <- NA
     d$up[nrow(d) - 1] <- NA
     d$lo[nrow(d) - 1] <- NA
     d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
     xint <- 0

     
     d <- d %>% mutate(outlier = ifelse(SNP %in% outliers_list, T, F))
     
     out <- ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
       ggplot2::geom_vline(xintercept = xint, linetype = "dotted") + 
       ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                            xmax = up, size = as.factor(outlier), colour = as.factor(outlier)), 
                               height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(outlier))) + 
       ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                             "")), colour = "grey") + 
       ggplot2::scale_colour_manual(values = c('FALSE' = "black",  'TRUE' = outliers_colour)) +
       ggplot2::scale_size_manual(values = c(0.3,  1)) + 
       ggplot2::theme_bw()+
       ggplot2::theme(legend.position = "none", 
                      axis.text.y = ggplot2::element_text(size = 8), 
                      axis.ticks.y = ggplot2::element_line(size = 0), 
                      axis.title.x = ggplot2::element_text(size = 8)) + 
       ggplot2::labs(y = "", 
                     subtitle = paste0("Single SNP forest plot with ",method," outliers"),
                     x = paste0("MR effect size for\n'", 
                                        d$exposure[1], "' on '", d$outcome[1], "'"))
     return(out)
   
  
}

mr_forest_plot_clusters <- function (d,  outliers_df, outliers_colour_list) {
  # this is a modified version of 2SMR function 
  
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  
  #if (sum(!grepl("All", d$SNP)) < 2) {
  # return(blank_plot("Insufficient number of SNPs"))}
  
  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
  am <- grep("All|cluster", d$SNP, value = TRUE)
  d$up <- d$b + 1.96 * d$se
  d$lo <- d$b - 1.96 * d$se
  d$tot <- 0.01
  d$tot[d$SNP %in% am] <- 1
  d$SNP <- as.character(d$SNP)
  nom <- d$SNP[!d$SNP %in% am]
  nom <- nom[order(d$b)]
  d <- rbind(d, d[nrow(d), ])
  d$SNP[nrow(d) - 1] <- ""
  d$b[nrow(d) - 1] <- NA
  d$up[nrow(d) - 1] <- NA
  d$lo[nrow(d) - 1] <- NA
  
  d <- d %>% left_join(outliers_df, by= c("SNP"="rsID")) 
  
  d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
  xint <- 0
  

  
  out <- ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
    ggplot2::geom_vline(xintercept = xint, linetype = "dotted") + 
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, colour = as.factor(cluster)), height = 0) +
    ggplot2::geom_point(ggplot2::aes(colour = as.factor(cluster))) + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour = "black") + 
    ggplot2::scale_colour_manual(values = outliers_colour_list) +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "none", 
                   axis.text.y = ggplot2::element_text(size = 8), 
                   axis.ticks.y = ggplot2::element_line(size = 0), 
                   axis.title.x = ggplot2::element_text(size = 8)) + 
    ggplot2::labs(y = "", 
                  subtitle = "Single SNP forest plot with MR-Clust clusters",
                  x = paste0("MR effect size for\n'", 
                                     d$exposure[1], "' on '", d$outcome[1], "'"))
  return(out)
  
  
}


calc_steiger <- function(harmonised, exposure_ss, outcome_ss, outcome_ncase = NA, outcome_ncontrol =NA){
  # assumed both traits are continuous, not binary
  
  harmonised$samplesize.exposure <- exposure_ss
  harmonised$samplesize.outcome <- outcome_ss 
  
  if (NA %in% harmonised$eaf.outcome){
    # if eaf is not available for outcome GWAS, we can't apply r_func below, so we will treat this data as continuous
    print("EAF not available; can't analyses outcome as binary for steiger filtering")
    outcome_ncase = NA
    outcome_ncontrol =NA
  }
  
  if (!is.na(outcome_ncase) & !is.na(outcome_ncontrol)){
    # if outcome case/control is provided, calculate prevalence and add everything to harmonised 
    # this is done for binary outcome analyses
    prevelance = outcome_ncase / (outcome_ncase + outcome_ncontrol)
    harmonised$ncase.outcome = outcome_ncase
    harmonised$ncontrol.outcome = outcome_ncontrol
    harmonised$prevelance.outcome = prevelance
    
    ## estimate r for binary or continuous exposures/outcomes
    harmonised <- r_func(harmonised, logistic.exposure=F, logistic.outcome=T)
  }
    
  
  harmonised <- steiger_filtering(harmonised) 
  harmonised_sub <- harmonised %>%  select(SNP, exposure, beta.exposure, pval.exposure, outcome, beta.outcome, pval.outcome, rsq.exposure, rsq.outcome, steiger_dir, steiger_pval)
  
  
  directionality <- directionality_test(harmonised)
  
  N = unique(harmonised$samplesize.exposure) #sample size
  K = length(harmonised$SNP) #number of SNPs
  total_r2 <- sum(harmonised$rsq.exposure) 
  Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
  
  
  summary <- directionality
  summary$Fst <- Fstat
  summary$total_r2 <- total_r2
  summary <- summary %>% select(-c("id.outcome", "id.exposure")) 
  
  
  return(list(Fstat = Fstat,
              total_r2 = total_r2,
              directionality= directionality,
              summary = summary,
              single_rsq = harmonised_sub))
  
  
}

## estimate r for binary or continuous exposures/outcomes
r_func <- function(x, logistic.exposure,logistic.outcome ){
  # function from: https://github.com/sjfandrews/MR_ADPhenome/blob/b64d8821dbf1546090f47e0642cc8092592cddc8/workflow/scripts/mr_SteigerTest.R#L13
  x$r.exposure <- if(logistic.exposure == TRUE){
    x %>%
      mutate(r.exposure = get_r_from_lor(x$beta.exposure, x$eaf.exposure, x$ncase.exposure, x$ncontrol.exposure, x$prevelance.exposure)) %>% pull(r.exposure)
  } else if(logistic.exposure == FALSE){
    x %>% mutate(r.exposure = get_r_from_pn(x$pval.exposure, x$samplesize.exposure)) %>% pull(r.exposure)
  }
  
  x$r.outcome <- if(logistic.outcome == TRUE){
    x %>% mutate(r.outcome = get_r_from_lor(x$beta.outcome, x$eaf.outcome, x$ncase.outcome, x$ncontrol.outcome, x$prevelance.outcome)) %>% pull(r.outcome)
  } else if(logistic.outcome == FALSE){
    x %>% mutate(r.outcome = get_r_from_pn(x$pval.outcome, x$samplesize.outcome)) %>% pull(r.outcome)
  }
  x
}


# PheWAS
get_pheno_assoc <- function(snp, ao_eu){
  # ao_eu is made like this:
  #ao_eu <-TwoSampleMR::available_outcomes() %>% filter(population == "European")
  
  # query phenoscanner
  res <- phenoscanner::phenoscanner(snpquery=snp)
  phenosc <- res$results 
  
  if (nrow(phenosc) > 0){
    phenosc<- phenosc %>%  as.data.frame() %>% 
      filter(p <= 5e-8) %>% 
      filter(ancestry == "European") %>% 
      arrange(p) %>% pull(trait) %>% unique()
  } else {  phenosc<-c() } # in case no res
  
  # query opengwas
  phwnogw<- ieugwasr::phewas(variants=snp, pval=5e-8)
  if (nrow(phwnogw) > 0){
    phwnogw <- phwnogw %>%  as.data.frame()  %>% 
      filter(id %in% ao_eu$id) %>% 
      arrange(p) %>% pull(trait) %>% unique()
  } else {  phwnogw<-c() }  # in case no res
  
  # merge sources
  out <- sort(unique(c(phenosc, phwnogw)))
  
}




### eQTL
### 
### 

#library(httr)
#library(jsonlite)
#library(GenomicRanges)
#library(biomaRt)
#
eqtl_for_snps <- function(variants){
  eqtl_df <- data.frame()
  for (i in 1:length(variants)){
    eqtl_df<-bind_rows(eqtl_df, get_eQTL_for_snp(variants[i]))
    #print(paste0("done: ", i))
  }
  gene_names <- map_ensembl_to_genename(eqtl_df$gene_id)
  eqtl_df<-left_join(eqtl_df, gene_names, 
                     by = c('gene_id' = 'ensembl_gene_id')) %>% 
    dplyr::select('hgnc_symbol', everything()) %>% 
    arrange(desc(median_tpm))
}

get_eQTL_for_snp <- function(variant){
  
  request = httr::GET(url = "http://www.ebi.ac.uk/eqtl/api/associations", 
                      query = list(
                        variant_id = variant,
                        size = 1000,
                        p_upper = 5e-8)
  )
  #print(variant)
  stopifnot(request$status_code==200) # 200 is good
  
  response = httr::content(request, as = "text", encoding = "UTF-8")
  variant_assoc = jsonlite::fromJSON(response, flatten = TRUE)$`_embedded`$associations
  
  if (length(variant_assoc)!=0){
    for (i in 1:length(variant_assoc)) {
      variant_assoc[[i]]<-variant_assoc[[i]] %>% purrr::discard(is.null) # drop any null items
    }
    variant_assoc_df <- bind_rows(variant_assoc) %>% dplyr::select(rsid, gene_id, qtl_group, pvalue, everything()) %>% arrange(pvalue)
  }else{
    return(data.frame())
  }
  
  nextq = jsonlite::fromJSON(response, flatten = TRUE)$`_links`
  if( any(names(nextq)=="next")){
    stop("API limits 1000 results, but there might be more -- need to investigate")
  }
  return(variant_assoc_df)
}


map_ensembl_to_genename <- function(values){
  
  #retrieve all genes with their GRCh37 coordinates from biomart
  mart_grch37 = biomaRt::useEnsembl(biomart="ensembl",GRCh=37)
  mart_grch37 = biomaRt::useDataset("hsapiens_gene_ensembl", mart_grch37)
  
  # retrieve gene symbols using biomart (eQTL Catalog returns ensembl)
  mart_query = biomaRt::getBM(mart=mart_grch37,
                              attributes=c("ensembl_gene_id","hgnc_symbol"),
                              filters= c("ensembl_gene_id"), 
                              values=unique(values))
  return(mart_query)
}



## function to call enrich fro pathways
## 
enrich_dbs<-function(gene_list, dbs, adjpval_filter = 0.05){
  
  enriched <- enrichr(gene_list, dbs)
  # flatten list into a table; handle empty tables
  for (db_name in names(enriched)){
    if( dim(enriched[[db_name]])[1] > 0){
      enriched[[db_name]]$db <- db_name
    } else  {
      # if it's empty, delete it
      enriched[[db_name]] <- NULL
    }
  }
  enriched_df<-bind_rows(enriched)
  if (dim(enriched_df)[1] > 0){
    enriched_df<- enriched_df %>%
      filter(Adjusted.P.value < adjpval_filter) %>% 
      separate(Overlap, into = c("found_genes", "total_genes"), sep="/", remove = F)%>% 
      arrange(Odds.Ratio)
  } else{
    enriched_df<-data.frame()
  }
}



