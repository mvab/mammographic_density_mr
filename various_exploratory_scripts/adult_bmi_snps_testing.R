# prep fro adjusted GWAS


adult_body_size_exp <- read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tophits/adult_bmi_tophits.tsv")
dim(adult_body_size_exp) # 173 - we asked for those

adult_bmi_ukb_exp <- read_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tophits/adult_bmi_cont_ukb_tophits.tsv")
dim(adult_bmi_ukb_exp) # 132 - we need these (clumpled from all)

length(intersect(adult_body_size_exp$SNP, adult_bmi_ukb_exp$SNP)) # 30 # from all we have these



adult_bmi_ukb_gwas <- vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tidy/adult_bmi_cont_ukb_GWAS_tidy_outcome.txt.gz")
snp_sub <- adult_bmi_ukb_gwas %>% filter(SNP %in% adult_body_size_exp$SNP) 
dim(snp_sub)# 126 (all) # extract adultd body size SNP (we have) that are also in BMI GWAS

# extract instuments from those
adul_bmi_inst <- snp_sub %>% filter(pval.outcome < 5e-8) %>% 
  convert_outcome_to_exposure() %>% 
  clump_data(., clump_r2 = 0.001)


adul_bmi_inst %>% write_tsv("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tophits/adult_bmi_cont_ukb_tophits_from_body_size_hits.tsv")
dim(adul_bmi_inst) # dim 79 --- this is the set 


# double check those 79 SNPs in MD data
mediator <- vroom(paste0(data_path_gwas,'dense_area_GWAS_tidy_outcome.txt.gz'))
mediator %>%  filter(SNP %in% adult_bmi_ukb_exp$SNP) %>% dim() # 37 / 132

mediator %>%  filter(SNP %in% adul_bmi_inst$SNP) %>% dim() # 77 / 79




# check adult BMI effect on BC 

## full SNP set (aka 30)

out<-quick_mr(adult_bmi_ukb_exp)
out # 0.8382296 0.7603641 0.924068

## subset SNP set
out<-quick_mr(adul_bmi_inst) #(aka 79)
out # 0.8288922 0.7279833 0.9437885




## select random set of 30

set.seed(1)

test1 <- sample_n(adul_bmi_inst, 30)
out<-quick_mr(test1)
out # 0.7643165 0.5978957 0.9770597


test2 <- sample_n(adul_bmi_inst, 30)
out<-quick_mr(test2)
out # 0.7857656 0.6175783 0.9997559

                  
test3 <- sample_n(adul_bmi_inst, 30)
out<-quick_mr(test3)
out # 0.9393332 0.7947666 1.110196







quick_mr <- function(exp_dat){
  
  out <- extract_outcome_data(
    snps = exp_dat$SNP,
    outcome = "ieu-a-1126")
  
  harmonised<- harmonise_data(exposure_dat = exp_dat, 
                              outcome_dat = out)
  
  res_early <- TwoSampleMR::mr(harmonised, method_list = c('mr_ivw')) %>% 
    split_outcome() %>% 
    split_exposure() %>% 
    generate_odds_ratios() %>% 
    select(exposure, outcome, or,  or_lci95 , or_uci95, nsnp)
  
}



