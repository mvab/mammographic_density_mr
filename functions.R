
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


calc_steiger <- function(harmonised, exposure_ss, outcome_ss ){
  # assumed both traits are continuous, not binary
  
  harmonised$samplesize.exposure <- exposure_ss
  harmonised$samplesize.outcome <- outcome_ss 
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

