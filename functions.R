
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


calc_steiger <- function(harmonised, exposure_ss, outcome_ss, outcome_ncase = NA, outcome_ncontrol =NA){
  # assumed both traits are continuous, not binary
  
  harmonised$samplesize.exposure <- exposure_ss
  harmonised$samplesize.outcome <- outcome_ss 
  
  if (!is.na(outcome_ncase) & !is.na(outcome_ncontrol)){
    # if outcome case/control is provided, calculate prevalence and add everything to harmonised 
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
  out <- unique(c(phenosc, phwnogw))
  
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



