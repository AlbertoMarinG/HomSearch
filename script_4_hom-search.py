""" Script to analyze RAD51 signals from different experiments """
__author__ = "Alberto Marin"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Alberto Marin"
####################################### IMPORT FUNCTIONS #######################################
import src.newchip as nc
import src.hic as hic
import os
####################################### DEFINE SOME PATHS #######################################
# Other paths
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]
tads_rao = "lib/4DNFIBKY9EG9_human_boundaries_strong.bed"
file_BLISS = "lib/Best_BLISS_Sites_HEK293_MultitargetGG.csv"
alnpath = "lib/psearch_hg38_align.csv"
# Paths to data
datadir = "/mnt/e/ZZ_DATA_LE_PAPER/Processed_Data/"
wapl_dir = datadir + "WAPL-AID_Rad51/"
profs_wapl = wapl_dir + "0_Profiles/"
os.makedirs(profs_wapl) if not os.path.exists(profs_wapl) else None
sschip_dir = datadir + "Strand_Spec_ChIP/"
profs_sschp = sschip_dir + "Profiles/"
os.makedirs(profs_sschp) if not os.path.exists(profs_sschp) else None
combidir = datadir + "ChIPs_WT/Combined/"
profs_dir = combidir + "Profiles/"
os.makedirs(profs_dir) if not os.path.exists(profs_dir) else None
# Paths to merged RAD51 ChIP-Seq files
bam_r51_3h = combidir + "RAD51_AluGG_EP_merged.bam"
bam_r51_no = combidir + "RAD51_NoEP_merged.bam"
bam_r51_novus_3h = combidir + "RAD51_Novus_EP_BR12.bam"
bam_r51_novus_no = combidir + "RAD51_Novus_NoEP_BR12.bam"
# Paths to ssChIP
sschp_no1 = sschip_dir + "BR1/ssChIP_NoEP2_hg38_final.bam"
sschp_no2 = sschip_dir + "BR2/ssChIP_NoEP3_hg38_final.bam"
sschp_ep1 = sschip_dir + "BR1/ssChIP_AluGG2_hg38_final.bam"
sschp_ep2 = sschip_dir + "BR2/ssChIP_AluGG3_hg38_final.bam"
out_ss_no1 = sschip_dir + "Split_Reads/ssChIP_NoEP2_hg38"
out_ss_no2 = sschip_dir + "Split_Reads/ssChIP_NoEP3_hg38"
out_ss_ep1 = sschip_dir + "Split_Reads/ssChIP_AluGG2_hg38"
out_ss_ep2 = sschip_dir + "Split_Reads/ssChIP_AluGG3_hg38"
# Paths to WAPL-AID, RAD51 ChIP-Seq
# Rep-1
r51_wapl_noaux_noep_1 = wapl_dir + "BR1/Rad51_WAPL_NoAux_NoEP_hg38_final.bam"
r51_wapl_noaux_ep_1 = wapl_dir + "BR1/Rad51_WAPL_NoAux_AluGG_hg38_final.bam"
r51_wapl_aux_noep_1 = wapl_dir + "BR1/Rad51_WAPL_Aux_NoEP_hg38_final.bam"
r51_wapl_aux_ep_1 = wapl_dir + "BR1/Rad51_WAPL_Aux_AluGG_hg38_final.bam"
# Rep-2
r51_wapl_noaux_noep_2 = wapl_dir + "BR2/NoAux-NoEP_hg38_final.bam"
r51_wapl_noaux_ep_2 = wapl_dir + "BR2/NoAux-EP_hg38_final.bam"
r51_wapl_aux_noep_2 = wapl_dir + "BR2/Aux-NoEP_hg38_final.bam"
r51_wapl_aux_ep_2 = wapl_dir + "BR2/Aux-EP_hg38_final.bam"
# gRNA parameters
AluGG = "CCTGTAGTCCCAGCTACTGG"
################################## QUANTIFY RAD51 SPREAD ################################################
""" Get RAD51 widths """
nc.get_span_width(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, bam_r51_3h, bam_r51_no, profs_dir + "widths_rad51_ontgt")
nc.get_span_width(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, bam_r51_novus_3h, bam_r51_novus_no, profs_dir + "widths_rad51_novus_ontgt")
################################## RAD51 Enrichment Across TAD boundaries ###############################
""" Get RAD51 profiles around TADs/random sites
nc.peak_profile_wide(hic.gen_TAD_bounds_around_cuts(tads_rao, file_BLISS, hg38, 100, 700000, 200000, 150000), 
                     hg38, bam_r51_no, profs_dir + "RAD51_TADs_Comb_NoEP_Best100_200to700", span_rad=150000, res=10000, wind_rad=20000)
nc.peak_profile_wide(hic.gen_TAD_bounds_around_cuts(tads_rao, file_BLISS, hg38, 100, 700000, 200000, 150000), 
                     hg38, bam_r51_3h, profs_dir + "RAD51_TADs_Comb_AluGG_Best100_200to700", span_rad=150000, res=10000, wind_rad=20000)
nc.peak_profile_wide(hic.gen_rands_around_cuts(file_BLISS, hg38, 100, 4, 700000, 200000, 150000), 
                     hg38, bam_r51_no, profs_dir + "RAD51_TADs_Comb_NoEP_Rands_Best100_200to700", span_rad=150000, res=10000, wind_rad=20000)
nc.peak_profile_wide(hic.gen_rands_around_cuts(file_BLISS, hg38, 100, 4, 700000, 200000, 150000), 
                     hg38, bam_r51_3h, profs_dir + "RAD51_TADs_Comb_AluGG_Rands_Best100_200to700", span_rad=150000, res=10000, wind_rad=20000)"""
""" Quantify drop in RAD51 at TADs/random sites
nc.peak_profile_wide(hic.gen_TAD_bounds_around_cuts(tads_rao, file_BLISS, hg38, 100, 700000, 200000, 
                                                    150000), hg38, bam_r51_no, profs_dir + "RAD51_TADs_Comb_NoEP_Best100_200to700_coarse", span_rad=150000, res=50000, wind_rad=50000)
nc.peak_profile_wide(hic.gen_TAD_bounds_around_cuts(tads_rao, file_BLISS, hg38, 100, 700000, 200000, 
                                                    150000), hg38, bam_r51_3h, profs_dir + "RAD51_TADs_Comb_AluGG_Best100_200to700_coarse", span_rad=150000, res=50000, wind_rad=50000)
nc.peak_profile_wide(hic.gen_rands_around_cuts(file_BLISS, hg38, 100, 4, 700000, 200000, 
                                               150000), hg38, bam_r51_no, profs_dir + "RAD51_TADs_Comb_NoEP_Rands_Best100_200to700_coarse", span_rad=150000, res=50000, wind_rad=50000)
nc.peak_profile_wide(hic.gen_rands_around_cuts(file_BLISS, hg38, 100, 4, 700000, 200000, 
                                               150000), hg38, bam_r51_3h, profs_dir + "RAD51_TADs_Comb_AluGG_Rands_Best100_200to700_coarse", span_rad=150000, res=50000, wind_rad=50000)"""
################################## Analyze Strand-Specific RAD51 ChIP-Seq ###############################
"""
# Split strand-specific ChIP-Seq reads into forward/reverse
nc.split_reads_sschip(sschp_no1, out_ss_no1)
nc.split_reads_sschip(sschp_no2, out_ss_no2)
nc.split_reads_sschip(sschp_ep1, out_ss_ep1)
nc.split_reads_sschip(sschp_ep2, out_ss_ep2)
# Broad Strand-Specific RAD51 profiles
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, sschp_no1,
                    profs_sschp + "NoEP1_All_1_5Mb", span_rad=1500000, res=5000, wind_rad=10000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, sschp_ep1,
                    profs_sschp + "AluGG1_All_1_5Mb", span_rad=1500000, res=5000, wind_rad=10000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_no1 + "_sorted_for.bam",
                    profs_sschp + "NoEP1_For_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_no1 + "_sorted_rev.bam",
                    profs_sschp + "NoEP1_Rev_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_ep1 + "_sorted_for.bam",
                    profs_sschp + "AluGG1_For_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_ep1 + "_sorted_rev.bam",
                    profs_sschp + "AluGG1_Rev_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, sschp_no2,
                    profs_sschp + "NoEP2_All_1_5Mb", span_rad=1500000, res=5000, wind_rad=10000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, sschp_ep2,
                    profs_sschp + "AluGG2_All_1_5Mb", span_rad=1500000, res=5000, wind_rad=10000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_no2 + "_sorted_for.bam",
                    profs_sschp + "NoEP2_For_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_no2 + "_sorted_rev.bam",
                    profs_sschp + "NoEP2_Rev_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_ep2 + "_sorted_for.bam",
                    profs_sschp + "AluGG2_For_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, out_ss_ep2 + "_sorted_rev.bam",
                    profs_sschp + "AluGG2_Rev_1_5Mb_fine", span_rad=1500000, res=1000, wind_rad=1000, norm_type=None, alignpam=False)"""
"""
# Narrow Strand-Specific RAD51 profiles
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_no1 + "_sorted_for.bam",
                    profs_sschp + "NoEP1_For_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_no1 + "_sorted_rev.bam",
                    profs_sschp + "NoEP1_Rev_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_ep1 + "_sorted_for.bam",
                    profs_sschp + "AluGG1_For_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_ep1 + "_sorted_rev.bam",
                    profs_sschp + "AluGG1_Rev_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, sschp_no1,
                    profs_sschp + "NoEP1_All_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, sschp_ep1,
                    profs_sschp + "AluGG1_All_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_no2 + "_sorted_for.bam",
                    profs_sschp + "NoEP2_For_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_no2 + "_sorted_rev.bam",
                    profs_sschp + "NoEP2_Rev_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_ep2 + "_sorted_for.bam",
                    profs_sschp + "AluGG2_For_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, out_ss_ep2 + "_sorted_rev.bam",
                    profs_sschp + "AluGG2_Rev_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, sschp_no2,
                    profs_sschp + "NoEP2_All_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 15000, AluGG), hg38, sschp_ep2,
                    profs_sschp + "AluGG2_All_15kb", span_rad=15000, res=50, wind_rad=100, norm_type=None, alignpam=False)"""
######################################### WAPL-AID RAD51 ChIP-Seq ######################################
""" WAPL-AID, RAD51 ChIP-Seq profiles 
# WAPL-AID, Rep-1
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_noaux_noep_1, wapl_dir + "RAD51_NoAux_NoEP_1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_noaux_ep_1, wapl_dir + "RAD51_NoAux_EP_1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_aux_noep_1, wapl_dir + "RAD51_Aux_NoEP_1", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_aux_ep_1, profs_wapl + "RAD51_Aux_EP_1", span_rad=1500000, res=5000, wind_rad=10000)
# WAPL-AID, Rep-2
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_noaux_noep_2, wapl_dir + "RAD51_NoAux_NoEP_2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_noaux_ep_2, wapl_dir + "RAD51_NoAux_EP_2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_aux_noep_2, wapl_dir + "RAD51_Aux_NoEP_2", span_rad=1500000, res=5000, wind_rad=10000)
nc.peak_profile_wide(nc.target_gen(alnpath, hg38, 1500000, AluGG), hg38, r51_wapl_aux_ep_2, profs_wapl + "RAD51_Aux_EP_2", span_rad=1500000, res=5000, wind_rad=10000)"""