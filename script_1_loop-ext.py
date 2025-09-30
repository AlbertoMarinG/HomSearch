""" Script to analyze WT Hi-C data, RAD21 siRNA Hi-C and Cohesin ChIP-Seq WT """
__author__ = "Alberto Marin"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Alberto Marin"
####################################### IMPORT FUNCTIONS #######################################
import src.hic as hic
import pandas as pd
import src.newchip as nc
import os
####################################### DEFINE SOME PATHS #######################################
#### INPUTS' PATHS ####
# General paths
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]
chrom_path = "lib/hg38.sizes" # Path to chrom sizes
file_BLISS = "lib/Best_BLISS_Sites_HEK293_MultitargetGG.csv"
jcr_dir = "/mnt/d/SOFTWARE/juicer/CPU/common/"
alnpath = "lib/psearch_hg38_align.csv"
jcr_path = jcr_dir + "juicer_tools.1.9.9_jcuda.0.8.jar"
file_MRE = "lib/Best_MRE_Sites_HEK293_MultitargetGG.csv"
# Data paths
datadir = "/mnt/e/ZZ_DATA_LE_PAPER/Processed_Data/"
chipsdir = datadir + "ChIPs_WT/"
profsdir = chipsdir + "Cohesin_Profiles/"
os.makedirs(profsdir) if not os.path.exists(profsdir) else None
# Paths to data, Hi-C large
hic_large_no1 = datadir + "Large_HiC/Separate_Reps/a_NoEP/HiC_NoEP_B1_merged.hic"
hic_large_no2 = datadir + "Large_HiC/Separate_Reps/a_NoEP/HiC_NoEP_B2_merged.hic"
hic_large_ep1 = datadir + "Large_HiC/Separate_Reps/b_Cas9/HiC_AluGG_B1_merged.hic"
hic_large_ep2 = datadir + "Large_HiC/Separate_Reps/b_Cas9/HiC_AluGG_B2_merged.hic"
hic_all_no = datadir + "Large_HiC/Merged_Reps/HiC_NoEP_merge_all.hic"
hic_all_ep = datadir + "Large_HiC/Merged_Reps/HiC_AluGG_merge_all.hic"
# Paths to Hi-C data, RAD21 siRNA. Rep-1
neg_noep1 = datadir + "RAD21_HiC/BR1/neg_no_30.hic"
neg_ep1 = datadir + "RAD21_HiC/BR1/neg_ep_30.hic"
r21_noep1 = datadir + "RAD21_HiC/BR1/rad_no_30.hic"
r21_ep1 = datadir + "RAD21_HiC/BR1/rad_ep_30.hic"
# Paths to Hi-C data, RAD21 siRNA. Rep-2
neg_noep2 = datadir + "RAD21_HiC/BR2/neg_no_30.hic"
neg_ep2 = datadir + "RAD21_HiC/BR2/neg_ep_30.hic"
r21_noep2 = datadir + "RAD21_HiC/BR2/rad_no_30.hic"
r21_ep2 = datadir + "RAD21_HiC/BR2/rad_ep_30.hic"
# Paths to Hi-C data, RAD21-AID. Merged replicates
hic_noaux_noep = datadir + "RAD21-AID/1_Hi-C/HiC_NoEP_NoAux_merged.hic"
hic_noaux_ep = datadir + "RAD21-AID/1_Hi-C/HiC_EP_NoAux_merged.hic"
hic_aux_noep = datadir + "RAD21-AID/1_Hi-C/HiC_NoEP_Aux_merged.hic"
hic_aux_ep = datadir + "RAD21-AID/1_Hi-C/HiC_EP_Aux_merged.hic"
# Paths to cohesin ChIP-Seq data, Rep-1
chp_r21_noep1 = chipsdir + "BR1/Rad21_NoEP_hg38_final.bam"
chp_r21_ep1 = chipsdir + "BR1/Rad21_EP_hg38_final.bam"
chp_psmc_noep1 = chipsdir + "BR1/pSMC1_NoEP_hg38_final.bam"
chp_psmc_ep1 = chipsdir + "BR1/pSMC1_EP_hg38_final.bam"
chp_nip_noep1 = chipsdir + "BR1/HEK_NoEP_NIPBL_hg38_final.bam"
chp_nip_ep1 = chipsdir + "BR1/HEK_AluGG_NIPBL_hg38_final.bam"
# Paths to cohesin ChIP-Seq data, Rep-2
chp_r21_noep2 = chipsdir + "BR2/Rad21_NoEP_hg38_final.bam"
chp_r21_ep2 = chipsdir + "BR2/Rad21_EP_hg38_final.bam"
chp_psmc_noep2 = chipsdir + "BR2/pSMC1_NoEP_hg38_final.bam"
chp_psmc_ep2 = chipsdir + "BR2/pSMC1_EP_hg38_final.bam"
chp_nip_noep2 = chipsdir + "BR2/NIPBL_Prot-Tech_NoEP_hg38_final.bam"
chp_nip_ep2 = chipsdir + "BR2/NIPBL_Prot-Tech_EP_hg38_final.bam"
# gRNA and loop extrusion parameters
AluGG = "CCTGTAGTCCCAGCTACTGG"
width = 1500000
res_LE = 50000
num_tgts = 100
######################################### GET DSB LOOP PLOTS ####################################
targets_hic = hic.read_tgts_peaks(file_MRE, hg38, res_LE, width, max_cter=100, mrg_tgts=0)
""" Large Hi-C """
# Rep-1
out_csv_rats_mre_lar1 = datadir + "Large_HiC/Separate_Reps/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-bliss_br1"
mat_wt1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, hic_large_no1, hic_large_ep1, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_wt1).to_csv(out_csv_rats_mre_lar1 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_wt1, out_csv_rats_mre_lar1, "EP 3h", "NoEP", (width/1e6))
# Rep-2
out_csv_rats_mre_lar2 = datadir + "Large_HiC/Separate_Reps/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_br2"
mat_wt2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, hic_large_no2, hic_large_ep2, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_wt2).to_csv(out_csv_rats_mre_lar2 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_wt1, out_csv_rats_mre_lar2, "EP 3h", "NoEP", (width/1e6))
""" Get 4C-like plots """
targets_gen = hic.read_tgts_peaks(file_MRE, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0)
[all_3h_lft, all_3h_rgt] = hic.cut_centered_4c(targets_hic, jcr_path, datadir, hic_all_no, hic_all_ep, res_LE, width, datadir + "Large_HiC/Merged_Reps/CutCent_4C_3h_100-mre", vp_bins=2, cut_bin=0)
""" Loop plots for RAD21 siRNA experiments """
# Neg, Rep-1
out_loops_neg_1 = datadir + "RAD21_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_neg"
mat_loops_neg_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, neg_noep1, neg_ep1, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_neg_1).to_csv(out_loops_neg_1 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_neg_1, out_loops_neg_1, "EP 3h", "NoEP Neg siRNA", (width/1e6))
# RAD21, Rep-1
out_loops_r21_1 = datadir + "RAD21_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_r21"
mat_loops_r21_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, r21_noep1, r21_ep1, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_r21_1).to_csv(out_loops_r21_1 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_r21_1, out_loops_r21_1, "EP 3h", "NoEP RAD21 siRNA", (width/1e6))
# Neg, Rep-2
out_loops_neg_2 = datadir + "RAD21_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_neg"
mat_loops_neg_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, neg_noep2, neg_ep2, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_neg_2).to_csv(out_loops_neg_2 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_neg_2, out_loops_neg_2, "EP 3h", "NoEP Neg siRNA", (width/1e6))
# RAD21, Rep-2
out_loops_r21_2 = datadir + "RAD21_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_r21"
mat_loops_r21_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, r21_noep2, r21_ep2, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_r21_2).to_csv(out_loops_r21_2 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_r21_2, out_loops_r21_2, "EP 3h", "NoEP RAD21 siRNA", (width/1e6))
""" Loop plots for RAD21-AID experiments """
# No-Auxin
out_loops_noaux = datadir + "RAD21-AID/1_Hi-C/out_mtrx_NoEPvsCas9_noaux"
mat_loops_noaux = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, hic_noaux_noep, hic_noaux_ep, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_noaux).to_csv(out_loops_noaux + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_noaux, out_loops_noaux, "EP 3h", "RAD21 NoAuxin", (width/1e6))
# Auxin
out_loops_aux = datadir + "RAD21-AID/1_Hi-C/out_mtrx_NoEPvsCas9_aux"
mat_loops_aux = hic.get_mean_HiC_ratio(targets_hic, jcr_path, datadir, hic_aux_noep, hic_aux_ep, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_loops_aux).to_csv(out_loops_aux + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_loops_aux, out_loops_aux, "EP 3h", "RAD21 Auxin", (width/1e6))
""" Plot cohesin profiles """
# Biological replicate-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_r21_noep1, profsdir + "RAD21_1_5Mb_NoEP_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_r21_ep1, profsdir + "RAD21_1_5Mb_AluGG_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_psmc_noep1, profsdir + "pSMC1_1_5Mb_NoEP_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_psmc_ep1, profsdir + "pSMC1_1_5Mb_AluGG_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_nip_noep1, profsdir + "NIPBL_1_5Mb_NoEP_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_nip_ep1, profsdir + "NIPBL_1_5Mb_AluGG_BR1", span_rad=1500000, res=5000, wind_rad=10000)
# Biological replicate-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_r21_noep2, profsdir + "RAD21_1_5Mb_NoEP_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_r21_ep2, profsdir + "RAD21_1_5Mb_AluGG_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_psmc_noep2, profsdir + "pSMC1_1_5Mb_NoEP_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_psmc_ep2, profsdir + "pSMC1_1_5Mb_AluGG_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_nip_noep2, profsdir + "NIPBL_1_5Mb_NoEP_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, chp_nip_ep2, profsdir + "NIPBL_1_5Mb_EP_BR2", span_rad=1500000, res=5000, wind_rad=10000)