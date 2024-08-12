library(lassosum)
library(data.table)
params<-commandArgs(TRUE)
cc_trait<-params[1]
prog_trait<-params[2]
ref_panel_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/refpanel_dir/eur_refpanel_hm3_all")
cc_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/",cc_trait,"_cc_eur_all_std.txt")
cc_sst_df<-fread(cc_sst_file)
file_expect1<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/cc_",cc_trait,"_lassosum_res.rds")
res_lassosum_cc<-lassosum.pipeline(cor = cc_sst_df$beta,chr = cc_sst_df$chr,pos = cc_sst_df$bpos,
                                   A1 = cc_sst_df$a1,A2 = cc_sst_df$a2,snp = cc_sst_df$snpid,
                                   ref.bfile = ref_panel_file,LDblocks = "EUR.hg19")
saveRDS(res_lassosum_cc,file_expect1)

prog_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/",prog_trait,"_prog_hm3_eur_train_std.txt")
prog_sst_df<-fread(prog_sst_file)
file_expect2<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/prog_",prog_trait,"_lassosum_res.rds")
res_lassosum_prog<-lassosum.pipeline(cor = prog_sst_df$beta,chr = prog_sst_df$chr,pos = prog_sst_df$bpos,
                                   A1 = prog_sst_df$a1,A2 = prog_sst_df$a2,snp = prog_sst_df$snpid,
                                   ref.bfile = ref_panel_file,LDblocks = "EUR.hg19")
saveRDS(res_lassosum_prog,file_expect2)

mtag_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/mtag_processed/mtag_",
                      cc_trait,"-",prog_trait,".txt")
mtag_sst_df<-fread(mtag_sst_file)
file_expect3<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/mtag_",prog_trait,"_lassosum_res.rds")
res_lassosum_mtag<-lassosum.pipeline(cor = mtag_sst_df$beta,chr = mtag_sst_df$chr,pos = mtag_sst_df$bpos,
                                     A1 = mtag_sst_df$a1,A2 = mtag_sst_df$a2,snp = mtag_sst_df$snpid,
                                     ref.bfile = ref_panel_file,LDblocks = "EUR.hg19")
saveRDS(res_lassosum_mtag,file_expect3)