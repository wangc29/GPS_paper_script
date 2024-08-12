devtools::load_all("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Rscript/TLPRS_source/TLPRS-main/")
library(data.table)
###
params<-commandArgs(TRUE)
cc_trait<-params[1]
prog_trait<-params[2]
base_method<-params[3]
###
valid_ped_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/val_Data_correct/",
                       prog_trait,"_tlprs_val.ped")
ref_panel_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/simulation_genomewide/genotype_by_chrom/refpanel_dir/eur_refpanel_hm3_all")
valid_genotype_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/val_Data_correct/val.",prog_trait)
prior_beta_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/tlprs_input/",
                        cc_trait,"_",base_method,"_prior_beta.txt")
train_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/tlprs_input/",
                       prog_trait,"_tlprs_train_sst.txt")

outfile_tlprs<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/tlprs_beta/",
                      prog_trait,"_tlprs_",base_method,"_res")

res_tlprs<-TL_PRS(ped_file = valid_ped_file,
                  Covar_name = "SEX",
                  Y_name = "PHENO",
                  Ytype = "B" ,
                  train_file = ref_panel_file,
                  test_file = valid_genotype_file,
                  sum_stats_file = prior_beta_file,
                  target_sumstats_file = train_sst_file,
                  LDblocks = "EUR.hg19",
                  outfile =outfile_tlprs ,cluster = NULL)
res_out_path<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/tlprs_beta/",
                     prog_trait,"_tlprs_",base_method,"_res_final.rds")
saveRDS(res_tlprs,res_out_path)