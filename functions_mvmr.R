require("purrr")
require("tidyr")
require("tibble")
require("TwoSampleMR")

get_mv_exposures <- function(exposure_list, full_gwas_list, clump_exposures=FALSE) {
  
  # Here we are using a re-written source code of `mv_extract_exposures` function from 2SMR package.
  # It was neccesary to do it, as it only works with MRBase database input at the moment (but we have external exposures)
  
  
  # Get effects of each instrument from each exposure
  # all tophit SNPs in both exposures, (clumped optionally)
  exposures <- exposure_list %>% purrr::reduce(bind_rows) 

  if (clump_exposures) {
    # ***optional*** : clump exposures 
    temp <- exposures
    temp$id.exposure <- 1
    temp <- clump_data(temp)
    #temp <- clump_data_local(temp)
    exposures <- filter(exposures, SNP %in% temp$SNP)
  }
  
  
  # merge exposures (in 'outcomes' ~ full gwas format)
  # extract all instruments from full GWAS of exposures
  for (i in 1:length(full_gwas_list)){
    full_gwas_list[[i]] <- full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  d1 <- full_gwas_list %>%
    purrr::reduce(bind_rows) %>% 
    distinct()

  # get ids
  id_exposure <- unique(d1$id.outcome) 
  
  # convert first trait to exposure format  -- exp1 is exposure
  tmp_exposure <- d1 %>% filter(id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()
  # keep other traits as outcome -- exp2+ are outcomes
  tmp_outcome <- d1 %>% filter(id.outcome != id_exposure[1])
  
  # Harmonise against the first trait
  d <- harmonise_data(exposure_dat = tmp_exposure, 
                      outcome_dat = tmp_outcome, action=2)
  
  # Only keep SNPs that are present in all
  snps_not_in_all <- d %>% dplyr::count(SNP)  %>% 
                    filter(n < length(exposure_list)-1) %>%
                    pull(SNP)
  d <- filter(d, !SNP %in% snps_not_in_all)

  
  # Subset and concat data
  
  # for exp1 get exposure cols
  dh1x <- d %>% filter(id.outcome == id.outcome[1]) %>% 
    select(SNP, contains("exposure"))
  # for exp2 get outcome cols
  dh2x <-d %>%  select(SNP, contains("outcome"))
  # rename outcome to exposure in these
  names(dh2x) <- gsub("outcome", "exposure", names(dh2x) )
  # join together (drop not needed cols)
  exposure_dat <- bind_rows(dh1x, dh2x) %>%  
    select(-c("samplesize.exposure" ,"mr_keep.exposure", "pval_origin.exposure")) %>% 
    distinct()
  
  return(exposure_dat)
}



## NB: this is a shortcut from exposure_dat in 2SMR to exposures_joined (from manual steps)
# `exposure_dat` is from  `exposure_dat <- rbind.fill(dh1x, dh2x) ...` above, 
# it is used to generate equivant of `exposures_joined` from manual steps, but here
# it is called `exposures_joined_auto`


#  function to convert 2SMR format into MVMR format
make_mvmr_input <- function(exposure_dat, outcome.id.mrbase=NULL, outcome.data=NULL) {
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (!is.null(outcome.id.mrbase )) {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (!is.null(outcome.data)){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}




tidy_mvmr_output <- function(mvmr_res) {
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate) %>% 
    rename(se="Std. Error") %>% 
    rename(pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}



clump_data_local <- function(dat, clump_r2 = 0.001, path = "/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/"){
  #https://github.com/MRCIEU/TwoSampleMR/issues/173  
  dat %>% 
    rename(rsid = SNP, 
           pval = pval.exposure,
           id = id.exposure) %>% 
    ieugwasr::ld_clump(
      dat = .,
      clump_r2 = clump_r2,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = paste0(path, "01_Data/reference/1kg.v3/EUR")) %>% 
    rename(SNP = rsid, 
           pval.exposure = pval,
           id.exposure = id) 
  
}


