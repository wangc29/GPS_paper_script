library(data.table)
cc_trait_vec<-c("SLE","RA")
prog_trait_vec<-c("SLE_ANA","RA_RF")
trait_pair_df<-cbind(cc_trait_vec,prog_trait_vec)
for(jj in 1:nrow(trait_pair_df)){
  cc_trait<-as.character(trait_pair_df[jj,1])
  prog_trait<-as.character(trait_pair_df[jj,2])
  cc_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/mtag_input/",cc_trait,"_cc_eur_all_std.txt")
  prog_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/mtag_input/",prog_trait,"_prog_eur_train_std.txt")
  mtag_out_prefix<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/mtag_output/mtag_",cc_trait,"-",prog_trait)
  
  runCMD<-paste0("mtag.py",
                 " --sumstats ",paste0(cc_sst_file,",",prog_sst_file),
                 " --out ",mtag_out_prefix,
                 " --n_min 0.0",
                 " --stream_stdout",
                 " --std_betas",
                 " --force")
  system(runCMD)
}
