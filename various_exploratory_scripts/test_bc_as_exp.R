clust1<-c("rs10995181","rs3821098","rs75197674")

clust2<-c("rs11040963", "rs11877925", "rs149689338", "rs1503613", "rs17196752", "rs1892368", "rs2642278", "rs335160", "rs3819405", "rs4897108", "rs6703250", "rs6715731", "rs67901221", "rs6885843", "rs8098548")

clust3<-c("rs11205303","rs11684853","rs17625845","rs6001984")



clust1_bc_exp <- extract_outcome_data(clust1,"ieu-a-1126") %>% convert_outcome_to_exposure()
clust2_bc_exp <- extract_outcome_data(clust2,"ieu-a-1126") %>% convert_outcome_to_exposure()
clust3_bc_exp <- extract_outcome_data(clust3,"ieu-a-1126") %>% convert_outcome_to_exposure()


outcome <- vroom("/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project4/01_Data/GWAS_tidy/dense_area_GWAS_tidy_outcome.txt.gz")


harm1 <- harmonise_data(clust1_bc_exp, outcome)
harm2 <- harmonise_data(clust2_bc_exp, outcome)
harm3 <- harmonise_data(clust3_bc_exp, outcome)


res1 <- mr(harm1) %>% 
  split_exposure() %>% 
  generate_odds_ratios() %>% filter(method =="Inverse variance weighted") %>% select(exposure, outcome,nsnp,pval, starts_with("or"))
res2 <- mr(harm2) %>% 
  split_exposure() %>% 
  generate_odds_ratios()%>% filter(method =="Inverse variance weighted")%>% select(exposure, outcome, nsnp,pval, starts_with("or"))
res3 <- mr(harm3) %>% 
  split_exposure() %>% 
  generate_odds_ratios()%>% filter(method =="Inverse variance weighted")%>% select(exposure, outcome,nsnp, pval, starts_with("or"))


rbind(res1,res2,res3)
oik