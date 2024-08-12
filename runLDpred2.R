library(data.table)
library(bigsnpr)
params<-commandArgs(TRUE)
cc_trait<-params[1]
prog_trait<-params[2]
##LD info
ld_mt<-readRDS("/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/1kg_resource/ldblk_1kg_eur/all_ld_1kg_eur_corrected_matrix.rds")
ld_score<-Matrix::colSums(ld_mt^2)
hm3_info<-readRDS("/gpfs/group/dxl46/default/private/chenwang/project/GPS/dataset/ukb_resource/LDpred2_LDrefs/Hapmap3/map.rds")
hm3_info[,"marker_id"]<-paste0(hm3_info$chr,"_",hm3_info$pos,"_",hm3_info$a0,"_",hm3_info$a1,"_b37")
##
cc_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/",cc_trait,"_cc_eur_all_std.txt")
cc_sst_df<-fread(cc_sst_file)
match_ix_cc<-match(cc_sst_df$marker_id,hm3_info$marker_id)
ldsc_res_cc <- snp_ldsc(ld_score =ld_score[match_ix_cc], length(ld_score[match_ix_cc]),
                        chi2 = (cc_sst_df[,"beta"]/cc_sst_df[,"beta_se"])^2,
                        sample_size = cc_sst_df[,3][1], blocks = NULL)
h2_seq_cc <- round(ldsc_res_cc[["h2"]] * c(0.3, 0.7, 1, 1.4), 4)
p_seq_cc <- signif(seq_log(1e-3, 0.1, length.out = 10), 2)
ldpred_params_cc <- expand.grid(p = p_seq_cc, h2 = h2_seq_cc, sparse = c(FALSE))
ldpred2_res_cc<-snp_ldpred2_grid(corr = as_SFBM(ld_mt[match_ix_cc,match_ix_cc]),df_beta = cc_sst_df,grid_param =ldpred_params_cc,ncores = 1)
##
prog_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/REAL_DATA/Data/processed_sst/prog_trait/train_data/",prog_trait,"_prog_hm3_eur_train_std.txt")
prog_sst_df<-fread(prog_sst_file)
match_ix_prog<-match(prog_sst_df$marker_id,hm3_info$marker_id)
ldsc_res_prog <- snp_ldsc(ld_score =ld_score[match_ix_prog], length(ld_score[match_ix_prog]),
                        chi2 = (prog_sst_df[,"beta"]/prog_sst_df[,"beta_se"])^2,
                        sample_size = prog_sst_df[,3][1], blocks = NULL)
h2_seq_prog <- round(ldsc_res_prog[["h2"]] * c(0.3, 0.7, 1, 1.4), 4)
p_seq_prog <- signif(seq_log(1e-3, 0.1, length.out = 10), 2)
ldpred_params_prog <- expand.grid(p = p_seq_prog, h2 = h2_seq_prog, sparse = c(FALSE))
ldpred2_res_prog<-snp_ldpred2_grid(corr = as_SFBM(ld_mt[match_ix_prog,match_ix_prog]),df_beta = prog_sst_df,grid_param =ldpred_params_prog,ncores = 1)
##
mtag_sst_file<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/mtag_processed/mtag_",
                      cc_trait,"-",prog_trait,".txt")
mtag_sst_df<-fread(mtag_sst_file)
match_ix_mtag<-match(mtag_sst_df$marker_id,hm3_info$marker_id)
ldsc_res_mtag <- snp_ldsc(ld_score =ld_score[match_ix_mtag], length(ld_score[match_ix_mtag]),
                        chi2 = (mtag_sst_df[,"beta"]/mtag_sst_df[,"beta_se"])^2,
                        sample_size = mtag_sst_df[,3][1], blocks = NULL)
h2_seq_mtag <- round(ldsc_res_mtag[["h2"]] * c(0.3, 0.7, 1, 1.4), 4)
p_seq_mtag <- signif(seq_log(1e-3, 0.1, length.out = 10), 2)
ldpred_params_mtag <- expand.grid(p = p_seq_mtag, h2 = h2_seq_mtag, sparse = c(FALSE))
ldpred2_res_mtag<-snp_ldpred2_grid(corr = as_SFBM(ld_mt[match_ix_mtag,match_ix_mtag]),df_beta = mtag_sst_df,grid_param =ldpred_params_mtag,ncores = 1)
##
outpath_ldpred_cc<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/cc_",cc_trait,"_ldpred_beta.rds")
outpath_ldpred_prog<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/prog_",prog_trait,"_ldpred_beta.rds")
outpath_ldpred_mtag<-paste0("/gpfs/group/dxl46/default/private/chenwang/project/GPS/Results/mtag_",prog_trait,"_ldpred_beta.rds")
saveRDS(ldpred2_res_cc,outpath_ldpred_cc)
saveRDS(ldpred2_res_prog,outpath_ldpred_prog)
saveRDS(ldpred2_res_mtag,outpath_ldpred_mtag)