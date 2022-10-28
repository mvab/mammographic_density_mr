
library(TwoSampleMR)
library(tidyverse)


extract_tophits <- function(outcome_gwas){
  outcome_gwas %>%
    filter(pval.outcome < 5e-8) %>% 
    convert_outcome_to_exposure() %>% 
    clump_data(., clump_r2 = 0.001)
}


BCAC_DA <- vroom("../01_Data/GWAS_tidy/BCAC_DA_GWAS_tidy_outcome.txt.gz")
gws_bcac <-BCAC_DA %>% filter(pval.outcome <= 5e-8) %>% select(SNP, beta.outcome, pval.outcome)
dim(gws_bcac)

Sieh_DA <- vroom("../01_Data/GWAS_tidy/dense_area_GWAS_tidy_outcome.txt.gz")
gws_sieh <-Sieh_DA %>% filter(pval.outcome <= 5e-8) %>% select(SNP, beta.outcome, pval.outcome)
dim(gws_sieh)

inst_sieh <- extract_tophits(Sieh_DA)
dim(inst_sieh)
inst_bcac <- extract_tophits(BCAC_DA)
dim(inst_bcac)

length(intersect(inst_sieh$SNP, inst_bcac$SNP))

da_new <- read_tsv(paste0(data_path_tophits, "dense_area_tophits.tsv")) # 25


gws_bcac %>% filter(SNP %in% Sieh_DA$SNP) %>% left_join(gws_sieh, by ="SNP") %>% 
  mutate(inst = ifelse(SNP %in%da_new$SNP, T,F) )%>%  View()



out_sieh <- extract_outcome_data(snps = inst_sieh$SNP, outcome = "ieu-a-1126")
harmonised_sieh<- harmonise_data(exposure_dat = inst_sieh, outcome_dat = out_sieh)

res_single_sieh <- mr_singlesnp(harmonised_sieh, all_method=c("mr_ivw_fe","mr_ivw_mre","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))%>% 
  mutate(SNP = gsub("Inverse variance weighted", "IVW", SNP))
p2_sieh <- mr_forest_plot(res_single_sieh)



out_bcac <- extract_outcome_data(snps = inst_bcac$SNP, outcome = "ieu-a-1126")
harmonised_bcac<- harmonise_data(exposure_dat = inst_bcac, outcome_dat = out_bcac)
res_single_bcac <- mr_singlesnp(harmonised_bcac, all_method=c("mr_ivw_fe","mr_ivw_mre","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))%>% 
  mutate(SNP = gsub("Inverse variance weighted", "IVW", SNP))
p2_bcac <- mr_forest_plot(res_single_bcac)


