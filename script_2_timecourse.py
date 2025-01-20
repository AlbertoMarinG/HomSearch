""" Script to analyze time-course experiments """
__author__ = "Alberto Marin"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Alberto Marin"
####################################### IMPORT FUNCTIONS #######################################
import src.hic as hic
import pandas as pd
import src.newchip as nc
####################################### DEFINE SOME PATHS #######################################
# Data paths
datadir = "/mnt/e/ZZ_DATA_LE_PAPER/Processed_Data/"
tcdir = datadir + "Time_Course/"
# Other paths
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]
chrom_path = "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/hg38.chrom.sizes" # Path to chrom sizes
file_BLISS = "lib/Best_BLISS_Sites_HEK293_MultitargetGG.csv"
jcr_dir = "/mnt/d/SOFTWARE/juicer/CPU/common/"
alnpath = "lib/psearch_hg38_align.csv"
jcr_path = jcr_dir + "juicer_tools.1.9.9_jcuda.0.8.jar"
file_MRE = "lib/Best_MRE_Sites_HEK293_MultitargetGG.csv"
# Paths to time-course Hi-C. hic matrices
# Rep-1
hic_tc1_0m = tcdir + "1_HiC/BR1/0m_inter_30.hic"
hic_tc1_15m = tcdir + "1_HiC/BR1/15m_inter_30.hic"
hic_tc1_30m = tcdir + "1_HiC/BR1/30m_inter_30.hic"
hic_tc1_1h = tcdir + "1_HiC/BR1/1h_inter_30.hic"
hic_tc1_3h = tcdir + "1_HiC/BR1/3h_inter_30.hic"
# Rep-2
hic_tc2_0m = tcdir + "1_HiC/BR2/0m_inter_30.hic"
hic_tc2_15m = tcdir + "1_HiC/BR2/15m_inter_30.hic"
hic_tc2_30m = tcdir + "1_HiC/BR2/30m_inter_30.hic"
hic_tc2_1h = tcdir + "1_HiC/BR2/1h_inter_30.hic"
hic_tc2_3h = tcdir + "1_HiC/BR2/3h_inter_30.hic"
# Paths to time-course Hi-C. Insulation scores
# Rep-1
is_tc1_0m = tcdir + "1_HiC/BR1/inscore_0m.dat"
is_tc1_15m = tcdir + "1_HiC/BR1/inscore_15m.dat"
is_tc1_30m = tcdir + "1_HiC/BR1/inscore_30m.dat"
is_tc1_1h = tcdir + "1_HiC/BR1/inscore_1h.dat"
is_tc1_3h = tcdir + "1_HiC/BR1/inscore_3h.dat"
# Rep-2
is_tc2_0m = tcdir + "1_HiC/BR2/inscore_0m.dat"
is_tc2_15m = tcdir + "1_HiC/BR2/inscore_15m.dat"
is_tc2_30m = tcdir + "1_HiC/BR2/inscore_30m.dat"
is_tc2_1h = tcdir + "1_HiC/BR2/inscore_1h.dat"
is_tc2_3h = tcdir + "1_HiC/BR2/inscore_3h.dat"
# Paths to time-course ChIPs
# yH2AX, Rep-1
gh2a_tc1_0m = tcdir + "2_gH2AX/BR1/gH2AX_0min_hg38_final.bam"
gh2a_tc1_15m = tcdir + "2_gH2AX/BR1/gH2AX_15min_hg38_final.bam"
gh2a_tc1_30m = tcdir + "2_gH2AX/BR1/gH2AX_30min_hg38_final.bam"
gh2a_tc1_1h = tcdir + "2_gH2AX/BR1/gH2AX_1h_hg38_final.bam"
gh2a_tc1_3h = tcdir + "2_gH2AX/BR1/gH2AX_3h_hg38_final.bam"
# yH2AX, Rep-2
gh2a_tc2_0m = tcdir + "2_gH2AX/BR2/gH2AX_0m_hg38_final.bam"
gh2a_tc2_15m = tcdir + "2_gH2AX/BR2/gH2AX_15m_hg38_final.bam"
gh2a_tc2_30m = tcdir + "2_gH2AX/BR2/gH2AX_30m_hg38_final.bam"
gh2a_tc2_1h = tcdir + "2_gH2AX/BR2/gH2AX_1h_hg38_final.bam"
gh2a_tc2_3h = tcdir + "2_gH2AX/BR2/gH2AX_3h_hg38_final.bam"
# 53BP1, Rep-1
bp53_tc1_0m = tcdir + "3_53BP1/BR1/no-light-53BP1_hg38_final.bam"
bp53_tc1_15m = tcdir + "3_53BP1/BR1/15m-light-53BP1_hg38_final.bam"
bp53_tc1_30m = tcdir + "3_53BP1/BR1/30m-light-53BP1_hg38_final.bam"
bp53_tc1_1h = tcdir + "3_53BP1/BR1/1h-light-53BP1_hg38_final.bam"
bp53_tc1_3h = tcdir + "3_53BP1/BR1/3h-light-53BP1_hg38_final.bam"
# 53BP1, Rep-2
bp53_tc2_0m = tcdir + "3_53BP1/BR2/53BP1_0m_hg38_final.bam"
bp53_tc2_15m = tcdir + "3_53BP1/BR2/53BP1_15m_hg38_final.bam"
bp53_tc2_30m = tcdir + "3_53BP1/BR2/53BP1_30m_hg38_final.bam"
bp53_tc2_1h = tcdir + "3_53BP1/BR2/53BP1_1h_hg38_final.bam"
bp53_tc2_3h = tcdir + "3_53BP1/BR2/53BP1_3h_hg38_final.bam"
# RPA, Rep-1
rpa_tc1_0m = tcdir + "4_RPA/BR1/RPA_0m_hg38_final.bam"
rpa_tc1_15m = tcdir + "4_RPA/BR1/RPA_15m_hg38_final.bam"
rpa_tc1_30m = tcdir + "4_RPA/BR1/RPA_30m_hg38_final.bam"
rpa_tc1_1h = tcdir + "4_RPA/BR1/RPA_1h_hg38_final.bam"
rpa_tc1_3h = tcdir + "4_RPA/BR1/RPA_3h_hg38_final.bam"
# RPA, Rep-2
rpa_tc2_0m = tcdir + "4_RPA/BR2/RPA_0m_hg38_final.bam"
rpa_tc2_15m = tcdir + "4_RPA/BR2/RPA_15m_hg38_final.bam"
rpa_tc2_30m = tcdir + "4_RPA/BR2/RPA_30m_hg38_final.bam"
rpa_tc2_1h = tcdir + "4_RPA/BR2/RPA_1h_hg38_final.bam"
rpa_tc2_3h = tcdir + "4_RPA/BR2/RPA_3h_hg38_final.bam"
# RAD51, Rep-1
rd51_tc1_0m = tcdir + "5_Rad51/BR1/Rad51_0m_hg38_final.bam"
rd51_tc1_15m = tcdir + "5_Rad51/BR1/Rad51_15m_hg38_final.bam"
rd51_tc1_30m = tcdir + "5_Rad51/BR1/Rad51_30m_hg38_final.bam"
rd51_tc1_1h = tcdir + "5_Rad51/BR1/Rad51_1h_hg38_final.bam"
rd51_tc1_3h = tcdir + "5_Rad51/BR1/Rad51_3h_hg38_final.bam"
# RAD51, Rep-2
rd51_tc2_0m = tcdir + "5_Rad51/BR2/Rad51_0m_hg38_final.bam"
rd51_tc2_15m = tcdir + "5_Rad51/BR2/Rad51_15m_hg38_final.bam"
rd51_tc2_30m = tcdir + "5_Rad51/BR2/Rad51_30m_hg38_final.bam"
rd51_tc2_1h = tcdir + "5_Rad51/BR2/Rad51_1h_hg38_final.bam"
rd51_tc2_3h = tcdir + "5_Rad51/BR2/Rad51_3h_hg38_final.bam"
#################################################################################################
######################################### GET DSB LOOP PLOTS ####################################
#################################################################################################
# gRNA and loop extrusion parameters
AluGG = "CCTGTAGTCCCAGCTACTGG"
width = 1500000
res_LE = 50000
num_tgts = 100
targets_hic = hic.read_tgts_peaks(file_MRE, hg38, res_LE, width, max_cter=100, mrg_tgts=0)
""" Hi-C tc-1 """
# 15 min
out_csv_rats_mre_15min = tcdir + "1_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_15min"
mat_15m_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc1_0m, hic_tc1_15m, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_15m_1).to_csv(out_csv_rats_mre_15min + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_15m_1, out_csv_rats_mre_15min, "15 min Light", "No Light", (width/1e6))
# 30 min
out_csv_rats_mre_30min = tcdir + "1_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_30min"
mat_30m_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc1_0m, hic_tc1_30m, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_30m_1).to_csv(out_csv_rats_mre_30min + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_30m_1, out_csv_rats_mre_30min, "30 min Light", "No Light", (width/1e6))
# 1 h
out_csv_rats_mre_1h = tcdir + "1_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_1h"
mat_1h_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc1_0m, hic_tc1_1h, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_1h_1).to_csv(out_csv_rats_mre_1h + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_1h_1, out_csv_rats_mre_1h, "1h Light", "No Light", (width/1e6))
# 3 h
out_csv_rats_mre_3h = tcdir + "1_HiC/BR1/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_3h"
mat_3h_1 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc1_0m, hic_tc1_3h, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_3h_1).to_csv(out_csv_rats_mre_3h + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_3h_1, out_csv_rats_mre_3h, "3h Light", "No Light", (width/1e6))
""" Hi-C tc-2 """
# 15 min
out_csv_rats_mre_15min = tcdir + "1_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_15min"
mat_15m_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc2_0m, hic_tc2_15m, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_15m_2).to_csv(out_csv_rats_mre_15min + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_15m_2, out_csv_rats_mre_15min, "15 min Light", "No Light", (width/1e6))
# 30 min
out_csv_rats_mre_30min = tcdir + "1_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_30min"
mat_30m_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc2_0m, hic_tc2_30m, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_30m_2).to_csv(out_csv_rats_mre_30min + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_30m_2, out_csv_rats_mre_30min, "30 min Light", "No Light", (width/1e6))
# 1 h
out_csv_rats_mre_1h = tcdir + "1_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_1h"
mat_1h_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc2_0m, hic_tc2_1h, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_1h_2).to_csv(out_csv_rats_mre_1h + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_1h_2, out_csv_rats_mre_1h, "1h Light", "No Light", (width/1e6))
# 3 h
out_csv_rats_mre_3h = tcdir + "1_HiC/BR2/out_mtrx_NoEPvsCas9_5kres_1_5Mbw_100-mre_3h"
mat_3h_2 = hic.get_mean_HiC_ratio(targets_hic, jcr_path, tcdir, hic_tc2_0m, hic_tc2_3h, res_LE, width, obs='oe', norm='KR')
pd.DataFrame(mat_3h_2).to_csv(out_csv_rats_mre_3h + ".csv", header=None, index=False)
hic.plot_LE_matrix(mat_3h_2, out_csv_rats_mre_3h, "3h Light", "No Light", (width/1e6))
#################################################################################################
################################## Time-course ChIP profiles ####################################
#################################################################################################
"""  """
# gH2AX, Rep-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc1_0m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_0m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc1_15m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_15m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc1_30m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_30m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc1_1h, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_1h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc1_3h, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_3h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
# gH2AX, Rep-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc2_0m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_0m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc2_15m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_15m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc2_30m, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_30m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc2_1h, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_1h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, gh2a_tc2_3h, tcdir + "2_gH2AX/" + "gH2AX_1_5Mb_3h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
# 53BP1, Rep-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc1_0m, tcdir + "3_53BP1/" + "BP53_1_5Mb_0m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc1_15m, tcdir + "3_53BP1/" + "BP53_1_5Mb_15m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc1_30m, tcdir + "3_53BP1/" + "BP53_1_5Mb_30m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc1_1h, tcdir + "3_53BP1/" + "BP53_1_5Mb_1h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc1_3h, tcdir + "3_53BP1/" + "BP53_1_5Mb_3h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
# 53BP1, Rep-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc2_0m, tcdir + "3_53BP1/" + "BP53_1_5Mb_0m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc2_15m, tcdir + "3_53BP1/" + "BP53_1_5Mb_15m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc2_30m, tcdir + "3_53BP1/" + "BP53_1_5Mb_30m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc2_1h, tcdir + "3_53BP1/" + "BP53_1_5Mb_1h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, bp53_tc2_3h, tcdir + "3_53BP1/" + "BP53_1_5Mb_3h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
# RPA, Rep-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc1_0m, tcdir + "4_RPA/" + "RPA_1_5Mb_0m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc1_15m, tcdir + "4_RPA/" + "RPA_1_5Mb_15m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc1_30m, tcdir + "4_RPA/" + "RPA_1_5Mb_30m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc1_1h, tcdir + "4_RPA/" + "RPA_1_5Mb_1h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc1_3h, tcdir + "4_RPA/" + "RPA_1_5Mb_3h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
# RPA, Rep-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc2_0m, tcdir + "4_RPA/" + "RPA_1_5Mb_0m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc2_15m, tcdir + "4_RPA/" + "RPA_1_5Mb_15m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc2_30m, tcdir + "4_RPA/" + "RPA_1_5Mb_30m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc2_1h, tcdir + "4_RPA/" + "RPA_1_5Mb_1h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rpa_tc2_3h, tcdir + "4_RPA/" + "RPA_1_5Mb_3h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
# Rad51, Rep-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc1_0m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_0m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc1_15m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_15m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc1_30m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_30m_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc1_1h, tcdir + "5_Rad51/" + "Rad51_1_5Mb_1h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc1_3h, tcdir + "5_Rad51/" + "Rad51_1_5Mb_3h_BR1", span_rad=1500000, res=5000, wind_rad=10000)
# Rad51, Rep-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc2_0m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_0m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc2_15m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_15m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc2_30m, tcdir + "5_Rad51/" + "Rad51_1_5Mb_30m_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc2_1h, tcdir + "5_Rad51/" + "Rad51_1_5Mb_1h_BR2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, rd51_tc2_3h, tcdir + "5_Rad51/" + "Rad51_1_5Mb_3h_BR2", span_rad=1500000, res=5000, wind_rad=10000)