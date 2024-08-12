library(data.table)
params<-commandArgs(TRUE)
cc_trait<-params[1]
prog_trait<-params[2]
chrom_num<-as.integer(params[3])
cc_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/prscs_input/",cc_trait,"_cc_eur_all_std.txt")
prog_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/prscs_input/",prog_trait,"_prog_eur_train_std.txt")
mtag_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/mtag_output/",
                      "prog_",prog_trait,"_mtag4prscs",".txt")
ld_ref_dir<-"/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/1kg_resource/ldblk_1kg_eur"
bim_prefix<-"/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/used_var"
cc_sample_size_df<-fread("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/cc_median_sample_size.txt")
prog_sample_size_df<-fread("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/prog_median_sample_size.txt")
N_cc<-cc_sample_size_df[cc_sample_size_df$Trait==cc_trait,2]
N_prog<-prog_sample_size_df[prog_sample_size_df$Trait==prog_trait,2]
####
for(phi_val in c(0.01,0.001,0.0001)){
  outpath_cc<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/PRS_CS_beta/prscs_beta_cc_",cc_trait)
  runCMD_cc<-paste0("PRScs.py",
                    " --ref_dir=",ld_ref_dir,
                    " --bim_prefix=",bim_prefix,
                    " --sst_file=",cc_sst_file,
                    " --n_gwas=",N_cc,
                    " --out_dir=",outpath_cc,
                    " --phi=",phi_val,
                    " --chrom=",chrom_num,
                    " --beta_std=","True")
  system(runCMD_cc)
  ####
  outpath_prog<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/PRS_CS_beta/prscs_beta_prog_",prog_trait)
  runCMD_prog<-paste0("PRScs.py",
                      " --ref_dir=",ld_ref_dir,
                      " --bim_prefix=",bim_prefix,
                      " --sst_file=",prog_sst_file,
                      " --n_gwas=",N_prog,
                      " --out_dir=",outpath_prog,
                      " --phi=",phi_val,
                      " --chrom=",chrom_num,
                      " --beta_std=","True")
  
  system(runCMD_prog)
  ####
  outpath_mtag<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/PRS_CS_beta/prscs_beta_mtag_",prog_trait)
  runCMD_mtag<-paste0("PRScs.py",
                      " --ref_dir=",ld_ref_dir,
                      " --bim_prefix=",bim_prefix,
                      " --sst_file=",mtag_sst_file,
                      " --n_gwas=",N_prog,
                      " --out_dir=",outpath_mtag,
                      " --phi=",phi_val,
                      " --chrom=",chrom_num,
                      " --beta_std=","True")
  system(runCMD_mtag)
}
