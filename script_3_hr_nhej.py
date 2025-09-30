""" Script to study the connection between loop extrusion and HR/NHEJ """
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
# Other paths
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]
chrom_path = "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/hg38.chrom.sizes" # Path to chrom sizes
file_BLISS = "lib/Best_BLISS_Sites_HEK293_MultitargetGG.csv"
jcr_dir = "/mnt/d/SOFTWARE/juicer/CPU/common/"
alnpath = "lib/psearch_hg38_align.csv"
jcr_path = jcr_dir + "juicer_tools.1.9.9_jcuda.0.8.jar"
file_MRE = "lib/Best_MRE_Sites_HEK293_MultitargetGG.csv"
# Data paths
datadir = "/mnt/e/ZZ_DATA_LE_PAPER/Processed_Data/"
chips_dir = datadir + "ChIPs_WT/"
profs_dir = chips_dir + "1_HR_Profiles/"
os.makedirs(profs_dir) if not os.path.exists(profs_dir) else None
# gRNA parameters
AluGG = "CCTGTAGTCCCAGCTACTGG"
# Paths to ChIP-Seq files. RAD51, RPA and Lig IV
rad51_no_1 = chips_dir + "BR1/Rad51_NoEP1_hg38_final.bam"
rad51_ep_1 = chips_dir + "BR1/Rad51_EP1_hg38_final.bam"
rad51_no_2 = chips_dir + "BR2/Rad51_NoEP2_hg38_final.bam"
rad51_ep_2 = chips_dir + "BR2/Rad51_EP2_hg38_final.bam"
rpa_no_1 = chips_dir + "BR1/RPA_NoEP1_hg38_final.bam"
rpa_ep_1 = chips_dir + "BR1/RPA_EP1_hg38_final.bam"
rpa_no_2 = chips_dir + "BR2/RPA_NoEP2_hg38_final.bam"
rpa_ep_2 = chips_dir + "BR2/RPA_EP2_hg38_final.bam"
rad51nov_no_1 = chips_dir + "BR1/Rad51Novus_NoEP_hg38_final.bam"
rad51nov_ep_1 = chips_dir + "BR1/Rad51Novus_EP1_hg38_final.bam"
rad51nov_no_2 = chips_dir + "BR2/HEK_NoEP_Rad51_Novus_hg38_final.bam"
rad51nov_ep_2 = chips_dir + "BR2/HEK_AluGG_Rad51_Novus_hg38_final.bam"
ligiv_no_1 = chips_dir + "BR1/LigIV_NoEP_hg38_final.bam"
ligiv_ep_1 = chips_dir + "BR1/LigIV_EP_hg38_final.bam"
ligiv_no_2 = chips_dir + "BR2/LigIV_NoEP2_hg38_final.bam"
ligiv_ep_2 = chips_dir + "BR2/LigIV_EP2_hg38_final.bam"
# Paths to ChIP-Seq files. RAD21, pSMC1 and NIPBL
chp_r21_noep1 = chips_dir + "BR1/Rad21_NoEP_hg38_final.bam"
chp_r21_ep1 = chips_dir + "BR1/Rad21_EP_hg38_final.bam"
chp_psmc_noep1 = chips_dir + "BR1/pSMC1_NoEP_hg38_final.bam"
chp_psmc_ep1 = chips_dir + "BR1/pSMC1_EP_hg38_final.bam"
chp_nip_noep1 = chips_dir + "BR1/HEK_NoEP_NIPBL_hg38_final.bam"
chp_nip_ep1 = chips_dir + "BR1/HEK_AluGG_NIPBL_hg38_final.bam"
chp_r21_noep2 = chips_dir + "BR2/Rad21_NoEP_hg38_final.bam"
chp_r21_ep2 = chips_dir + "BR2/Rad21_EP_hg38_final.bam"
chp_psmc_noep2 = chips_dir + "BR2/pSMC1_NoEP_hg38_final.bam"
chp_psmc_ep2 = chips_dir + "BR2/pSMC1_EP_hg38_final.bam"
chp_nip_noep2 = chips_dir + "BR2/NIPBL_Prot-Tech_NoEP_hg38_final.bam"
chp_nip_ep2 = chips_dir + "BR2/NIPBL_Prot-Tech_EP_hg38_final.bam"
# Paths to Hi-C cell cycle
hic_noli_g11 = datadir + "HiC_Cell_Cycle/1_G1/BR1/HiC_NoLight_30.hic"
hic_3hli_g11 = datadir + "HiC_Cell_Cycle/1_G1/BR1/HiC_3hLight_30.hic"
hic_noli_g21 = datadir + "HiC_Cell_Cycle/2_G2/BR1/HiC_NoLight_30.hic"
hic_3hli_g21 = datadir + "HiC_Cell_Cycle/2_G2/BR1/HiC_3hLight_30.hic"
hic_noli_g12 = datadir + "HiC_Cell_Cycle/1_G1/BR2/HiC_NoLight_30.hic"
hic_3hli_g12 = datadir + "HiC_Cell_Cycle/1_G1/BR2/HiC_3hLight_30.hic"
hic_noli_g22 = datadir + "HiC_Cell_Cycle/2_G2/BR2/HiC_NoLight.hic"
hic_3hli_g22 = datadir + "HiC_Cell_Cycle/2_G2/BR2/HiC_3hLight.hic"
# Paths to Hi-C drugs
hic_noep_mir1 = datadir + "HiC_Drugs/Mirin/BR1/inter_NoEP.hic"
hic_cas9_mir1 = datadir + "HiC_Drugs/Mirin/BR1/inter_AluGG.hic"
hic_noep_mir2 = datadir + "HiC_Drugs/Mirin/BR2/inter_NoEP_30.hic"
hic_cas9_mir2 = datadir + "HiC_Drugs/Mirin/BR2/inter_AluGG_30.hic"
# Paths to inscore, Mirin
is_mir_no1 = datadir + "HiC_Drugs/Mirin/BR1/inscore_noep.dat"
is_mir_ep1 = datadir + "HiC_Drugs/Mirin/BR1/inscore_ep.dat"
is_mir_no2 = datadir + "HiC_Drugs/Mirin/BR2/inscore_noep.dat"
is_mir_ep2 = datadir + "HiC_Drugs/Mirin/BR2/inscore_ep.dat"
# Paths to inscore outputs, Mirin
outis_mir_no_b1 = datadir + "HiC_Drugs/Mirin/BR1/On-tgt_InsScore_noep_b1"
outis_mir_ep_b1 = datadir + "HiC_Drugs/Mirin/BR1/On-tgt_InsScore_ep_b1"
outis_mir_no_b2 = datadir + "HiC_Drugs/Mirin/BR2/On-tgt_InsScore_noep_b2"
outis_mir_ep_b2 = datadir + "HiC_Drugs/Mirin/BR2/On-tgt_InsScore_ep_b2"
# Paths to inscore, large
is_no1 = datadir + "Large_HiC/Separate_Reps/a_NoEP/inscore_noep1.dat"
is_ep1 = datadir + "Large_HiC/Separate_Reps/b_Cas9/inscore_ep1.dat"
is_no2 = datadir + "Large_HiC/Separate_Reps/a_NoEP/inscore_noep2.dat"
is_ep2 = datadir + "Large_HiC/Separate_Reps/b_Cas9/inscore_ep2.dat"
# Paths to inscore outputs, large
outis_no_b1 = datadir + "Large_HiC/Separate_Reps/a_NoEP/On-tgt_InsScore_noep_b1"
outis_ep_b1 = datadir + "Large_HiC/Separate_Reps/b_Cas9/On-tgt_InsScore_ep_b1"
outis_no_b2 = datadir + "Large_HiC/Separate_Reps/a_NoEP/On-tgt_InsScore_noep_b2"
outis_ep_b2 = datadir + "Large_HiC/Separate_Reps/b_Cas9/On-tgt_InsScore_ep_b2"
outis_bls_no_b1 = datadir + "Large_HiC/Separate_Reps/a_NoEP/InsScore_noep_b1"
outis_bls_ep_b1 = datadir + "Large_HiC/Separate_Reps/b_Cas9/InsScore_ep_b1"
outis_bls_no_b2 = datadir + "Large_HiC/Separate_Reps/a_NoEP/InsScore_noep_b2"
outis_bls_ep_b2 = datadir + "Large_HiC/Separate_Reps/b_Cas9/InsScore_ep_b2"
# gRNA and loop extrusion parameters
width = 1500000
res_LE = 50000
""" Large, WT Hi-C
# Insulation score for each on-target site
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_no1, outis_no_b1)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_ep1, outis_ep_b1)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_no2, outis_no_b2)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_ep2, outis_ep_b2)
# Insulation score per cut-site, BLISS
hic.per_site_is_peaks(hic.read_tgts_peaks(file_BLISS, hg38, res=1, width=1000000, max_cter=100), is_no1, outis_bls_no_b1)
hic.per_site_is_peaks(hic.read_tgts_peaks(file_BLISS, hg38, res=1, width=1000000, max_cter=100), is_ep1, outis_bls_ep_b1)
hic.per_site_is_peaks(hic.read_tgts_peaks(file_BLISS, hg38, res=1, width=1000000, max_cter=100), is_no2, outis_bls_no_b2)
hic.per_site_is_peaks(hic.read_tgts_peaks(file_BLISS, hg38, res=1, width=1000000, max_cter=100), is_ep2, outis_bls_ep_b2)
# Insulation Score profiles
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_no1, outis_no_b1, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_ep1, outis_ep_b1, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_no2, outis_no_b2, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_ep2, outis_ep_b2, res = 25, span=20) """
""" Drugs Hi-C, Loop extrusion plots 
targets_gen = hic.read_tgts_peaks(file_MRE, hg38, res_LE, width, max_cter=100, mrg_tgts=0) """
"""
# Mirin, Rep-1
out_csv_rats_mre_mir1 = datadir + "HiC_Drugs/Mirin/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_mir1 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noep_mir1, hic_cas9_mir1, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_mir1).to_csv(out_csv_rats_mre_mir1 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_mir1, out_csv_rats_mre_mir1, "3h EP", "No EP", (width/1e6))
# Mirin, Rep-2
out_csv_rats_mre_mir2 = datadir + "HiC_Drugs/Mirin/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_mir2 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noep_mir2, hic_cas9_mir2, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_mir2).to_csv(out_csv_rats_mre_mir2 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_mir2, out_csv_rats_mre_mir2, "3h EP", "No EP", (width/1e6)) """
""" Drugs Hi-C, Insulation Scores
# Insulation Score profiles
# Mirin
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_no1, outis_mir_no_b1, span=1)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_ep1, outis_mir_ep_b1, span=1)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_no2, outis_mir_no_b2, span=1)
hic.per_site_is(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_ep2, outis_mir_ep_b2, span=1)
# Insulation Score profiles
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_no1, outis_mir_no_b1, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_ep1, outis_mir_ep_b1, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_no2, outis_mir_no_b2, res = 25, span=20)
hic.is_profiles(nc.target_gen(alnpath, hg38, 1, AluGG), is_mir_ep2, outis_mir_ep_b2, res = 25, span=20) """
""" Cell Cycle Hi-C
# G1, Rep-1
out_csv_rats_mre_g11 = datadir + "HiC_Cell_Cycle/1_G1/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_g1_1 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noli_g11, hic_3hli_g11, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_g1_1).to_csv(out_csv_rats_mre_g11 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_g1_1, out_csv_rats_mre_g11, "EP 3h-light", "No-light", (width/1e6))
# G2, Rep-1
out_csv_rats_mre_g21 = datadir + "HiC_Cell_Cycle/2_G2/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_g2_1 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noli_g21, hic_3hli_g21, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_g2_1).to_csv(out_csv_rats_mre_g21 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_g2_1, out_csv_rats_mre_g21, "EP 3h-light", "No-light", (width/1e6))
# G1, Rep-2
out_csv_rats_mre_g12 = datadir + "HiC_Cell_Cycle/1_G1/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_g1_2 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noli_g12, hic_3hli_g12, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_g1_2).to_csv(out_csv_rats_mre_g12 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_g1_2, out_csv_rats_mre_g12, "EP 3h-light", "No-light", (width/1e6))
# G2, Rep-2
out_csv_rats_mre_g22 = datadir + "HiC_Cell_Cycle/2_G2/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre"
mat_g2_2 = hic.get_mean_HiC_ratio(targets_gen, jcr_path, datadir, hic_noli_g22, hic_3hli_g22, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_g2_2).to_csv(out_csv_rats_mre_g22 + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_g2_2, out_csv_rats_mre_g22, "EP 3h-light", "No-light", (width/1e6)) """
""" ChIP-Seq profiles 
# Rad51, Reps 1, 2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_no_1, profs_dir + "Rad51_15kb_no_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_no_2, profs_dir + "Rad51_15kb_no_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_ep_1, profs_dir + "Rad51_15kb_ep_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_ep_2, profs_dir + "Rad51_15kb_ep_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_no_1, profs_dir + "Rad51_1_5Mb_no_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_no_2, profs_dir + "Rad51_1_5Mb_no_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_ep_1, profs_dir + "Rad51_1_5Mb_ep_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51_ep_2, profs_dir + "Rad51_1_5Mb_ep_BR2", span_rad=1500000, res=5000, wind_rad=10000)
# RPA, Reps 1, 2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rpa_no_1, profs_dir + "RPA_15kb_no_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rpa_no_2, profs_dir + "RPA_15kb_no_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rpa_ep_1, profs_dir + "RPA_15kb_ep_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rpa_ep_2, profs_dir + "RPA_15kb_ep_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_no_1, profs_dir + "RPA_1_5Mb_no_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_ep_1, profs_dir + "RPA_1_5Mb_ep_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_no_2, profs_dir + "RPA_1_5Mb_no_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_ep_2, profs_dir + "RPA_1_5Mb_ep_BR2", span_rad=1500000, res=5000, wind_rad=10000)"""
# Rad51 Novus, Reps 1, 2
""" 
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rad51nov_no_1, profs_dir + "Rad51novus_1_5Mb_no_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rad51nov_ep_1, profs_dir + "Rad51novus_1_5Mb_ep_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rad51nov_no_2, profs_dir + "Rad51novus_1_5Mb_no_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rad51nov_ep_2, profs_dir + "Rad51novus_1_5Mb_ep_BR2", span_rad=1500000, res=5000, wind_rad=10000)

nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51nov_no_1, profs_dir + "Rad51novus_15kb_no_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51nov_ep_1, profs_dir + "Rad51novus_15kb_ep_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51nov_no_2, profs_dir + "Rad51novus_15kb_no_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, rad51nov_ep_2, profs_dir + "Rad51novus_15kb_ep_BR2", span_rad=15000, res=50, wind_rad=100)"""
""" 
# DNA ligase IV, Reps 1, 2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, ligiv_no_1, profs_dir + "LigIV_15kb_no_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, ligiv_no_2, profs_dir + "LigIV_15kb_no_BR2", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, ligiv_ep_1, profs_dir + "LigIV_15kb_ep_BR1", span_rad=15000, res=50, wind_rad=100)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, ligiv_ep_2, profs_dir + "LigIV_15kb_ep_BR2", span_rad=15000, res=50, wind_rad=100)"""
""" Per-ste enrichment for correlation analysis
# Rep-1
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, rad51_ep_1, rad51_no_1, 10000, chips_dir + "BR1/Rad51_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, ligiv_ep_1, ligiv_no_1, 500, chips_dir + "BR1/LigIV_500bp_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_r21_ep1, chp_r21_noep1, 10000, chips_dir + "BR1/RAD21_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_psmc_ep1, chp_psmc_noep1, 10000, chips_dir + "BR1/pSMC1_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_nip_ep1, chp_nip_noep1, 10000, chips_dir + "BR1/NIPBL_10kb_100-bls")
# Rep-2
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, rad51_ep_2, rad51_no_2, 10000, chips_dir + "BR2/Rad51_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, ligiv_ep_2, ligiv_no_2, 500, chips_dir + "BR2/LigIV_500bp_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_r21_ep2, chp_r21_noep2, 10000, chips_dir + "BR2/RAD21_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_psmc_ep2, chp_psmc_noep2, 10000, chips_dir + "BR2/pSMC1_10kb_100-bls")
nc.per_site_enrchment(hic.read_tgts_peaks(file_BLISS, hg38, res_LE, width, max_cter=100, mrg_tgts=0, sort=1, cut_bin=0),
                      hg38, chp_nip_ep2, chp_nip_noep2, 10000, chips_dir + "BR2/NIPBL_10kb_100-bls") """