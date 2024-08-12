library(multivariateLassosum)
library(bigsnpr)
library(data.table)
params<-commandArgs(TRUE)
cc_trait<-params[1]
prog_trait<-params[2]
lambda_choice_num<-as.integer(params[3])
###
hm3_info<-readRDS("/storage/group/dxl46/default/private/chenwang/project/GPS/dataset/ukb_resource/LDpred2_LDrefs/Hapmap3/map.rds")
hm3_info[,"marker_id"]<-paste0(hm3_info$chr,"_",hm3_info$pos,"_",hm3_info$a0,"_",hm3_info$a1,"_b37")
cc_sst_file<-paste0("/storage/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/cc_trait/",cc_trait,"_cc_eur_all_std.txt")
cc_sst_df<-fread(cc_sst_file)
cc_sst_df[,"marker_id"]<-paste0(cc_sst_df$chr,"_",cc_sst_df$bpos,"_",cc_sst_df$a2,"_",cc_sst_df$a1,"_b37")
###
prog_sst_file<-paste0("/storage/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/",prog_trait,"_prog_hm3_eur_train_std.txt")
prog_sst_df<-fread(prog_sst_file)
prog_sst_df[,"marker_id"]<-paste0(prog_sst_df$chr,"_",prog_sst_df$bpos,"_",prog_sst_df$a2,"_",prog_sst_df$a1,"_b37")
##
cor_cc <- p2cor(p =cc_sst_df$pval, n = median(cc_sst_df$n,na.rm = T), sign = cc_sst_df$beta)
cor_prog <- p2cor(p =prog_sst_df$pval, n = median(prog_sst_df$n,na.rm = T), sign = prog_sst_df$beta)
###
###
Var.phenotypic <- c(1,1)
##Load genetic covariance
omega_file<-paste0("/storage/group/dxl46/default/private/chenwang/project/GPS/Data/sst_data/processed_sst/prog_trait/train_data/mtag_output/mtag_",
                   cc_trait,"-",prog_trait,"_omega_hat.txt")
phenotypic.genetic.Var.Cov.matrix <-as.matrix(read.table(omega_file))
##
ref_panel_file<-paste0("/storage/group/dxl46/default/private/chenwang/project/GPS/dataset/simulation_genomewide/genotype_by_chrom/refpanel_dir/eur_refpanel_hm3_all")
LDblocks <- "EUR.hg19"
lambda_vec_default<-exp(seq(log(0.001), log(1000), length.out=20))
outMulti <- multivariateLassosum::lassosum.pipeline(cor = list(cor_cc, cor_prog),
                                                    phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix,
                                                    Var.phenotypic = Var.phenotypic,
                                                    chr = list(cc_sst_df$chr, prog_sst_df$chr),
                                                    pos = list(cc_sst_df$bpos, prog_sst_df$bpos),
                                                    A1 = list(cc_sst_df$a1, prog_sst_df$a1),
                                                    A2 = list(cc_sst_df$a2, prog_sst_df$a2),
                                                    sample_size = c(median(cc_sst_df$n), median(prog_sst_df$n)),
                                                    ref.bfile = ref_panel_file,lambda = lambda_vec_default[lambda_choice_num],
                                                    LDblocks = LDblocks)

##Save the output
outpath<-paste0("/storage/group/dxl46/default/private/chenwang/project/GPS/Results/mvlassosum_beta/",
                prog_trait,"_mvl_raw_res_lambda_choice_",lambda_choice_num,".rds")
saveRDS(outMulti,outpath)