""" Helper script to analyze ChIP-Seq data to study repair of Cas9 breaks """
__author__ = "Alberto Marin"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Alberto Marin"
####################################### IMPORT FUNCTIONS #######################################
import pysam
import re
from collections import defaultdict
import numpy as np
import csv
import os

CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX', 'chrY']

def load_nparray(array):
    """ Load dataset as numpy array. """
    return np.loadtxt(array, dtype=object, delimiter=',')

def get_genome_dict(genome_str):
    """ Return dict that holds the number of base pairs for each chromosome in 'hg38','hg19','mm10'.

    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :return: dict with keys as chromosomes, values as maximum coordinate of each chromosome key
    """
    d = {}
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s.sizes" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            d[row[0]] = int(row[1])
    return d

def target_gen(alignfile, genome, span_r, guide):
    """ Generator to yield all putative on-target sites for a given protospacer

    :param alignfile: CSV file generated from get_targets_stats() output
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param guide: on-target protospacer sequence (no PAM)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM, discovered protospacer,
              # mismatches, non-mismatched protospacer )
            - The format is designed to be consistent with other generators that may yield
              mismatched sequences. Note that no mismatched sequences will be outputted by design.

    """
    aln = load_nparray(alignfile)
    g_dict = get_genome_dict(genome[0])
    pam_i = 'NGG'
    outlist = []
    for i in range(aln.shape[0]):
        row = aln[i, :]
        if row[0] == guide + pam_i:
            chr_i = row[4]
            if chr_i not in CHR:
                continue
            sen_i = row[6]
            if sen_i == '+':
                cut_i = int(row[5]) + 16
            else:
                cut_i = int(row[5]) + 6
            span_sta = max(1, cut_i - span_r)
            span_end = min(g_dict[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            mis_i = 0
            outlist.append([chr_i, span_rs, cut_i, sen_i, pam_i, guide, mis_i, guide])
    outlist = sorted(outlist, key=lambda x: (x[2]))
    outlist = sorted(outlist, key=lambda x: (x[0]))
    for out in outlist:
        yield tuple(out[1:])

def peak_profile_wide(generator, genome, bamfilein, fileout, span_rad=2000000, res=1000, wind_rad=10000,
                      norm_type=None, alignpam=True, wig=False):
    """ For each target location from generator, calculates enrichment at specified 'resolution'
        with sliding window of specified 'radius'. Outputs the enrichment from a BAM file as:
        (1) CSV file with each row one target location, column is enrichment values in a window
        centered at each target location, at a specific genomic resolution.
        (2) WIG file with the local (window-bounded) enrichment at each target location

    :param bamfilein: path to input BAM file
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param fileout: path to output file name (excludes extension)
    :param norm_type: If None (default), then normalize to RPM. If False, then no normalization
                       Otherwise, if list, then assume list of region strings for normalization,
                       i.e. ['chr12:6532000-6536000', 'chr15:44709500-44713500']
                       Otherwise, assume it is an alternative BAM file to use for RPM normalization.
    :param span_rad: radius of analysis, centered at the cut site | default 2E6 bp
    :param res: resolution, i.e. bp to skip to calculate enrichment | default 1E3 bp
    :param wind_rad: radius of sliding window | default 1E4 bp

    Results in two files (WIG and CSV), described above.
    """
    hgsize = get_genome_dict(genome[0])
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    chr_old, csv_peaks = None, []
    wlist_all = []
    numrows = int(span_rad * 2 / res) + 1
    wlist = None
    if norm_type is None:
        norm_num = bamin.mapped / 1E6
    elif not norm_type:
        norm_num = 1
    elif isinstance(norm_type, list):
        norm_num = sum(bamin.count(region=co) for co in norm_type) / len(norm_type) / 10
    else:
        bamalt = pysam.AlignmentFile(norm_type, 'rb')
        norm_num = bamalt.mapped / 1E6
        bamalt.close()
    for rs, cut, sen, pam, gui, mis, guide in generator:
        sense = int(sen+str(1)) if alignpam is True else 1
        chr_i = re.split('[:-]', rs)[0]
        sta_i = cut - span_rad * sense
        end_i = cut + span_rad * sense
        if min(sta_i, end_i) - wind_rad >= 0 and max(end_i, sta_i) + wind_rad < hgsize[chr_i]:
            wlist = [0] * numrows
            for row_i in range(numrows):
                center = sta_i + row_i * res * sense
                rs_i = "%s:%i-%i" % (chr_i, center - wind_rad, center + wind_rad)
                wlist[row_i] = bamin.count(region=rs_i) / norm_num
            wlist_all.append([chr_i, sta_i] + wlist)
            csv_peaks.append([chr_i, cut, gui + pam, mis] + wlist)
    bamin.close()
    # Compute the average enrichment and include it in the Excel file
    aver_list = []
    if wlist is not None:
        for col in range(len(wlist)):
            aver_i = 0
            for line in range(len(csv_peaks)):
                aver_i += float(csv_peaks[line][col + 4])
            aver_i = aver_i / int(len(csv_peaks))
            aver_list.append(aver_i)
    csv_peaks.append(["Average", "NA", "NA", "NA"] + aver_list)
    #
    if wig == True:
        _peak_profile_helper(wlist_all, res, fileout)
    head = ",".join(["chr", "cut", "guide+PAM", "mismatches"] +
                    list(map(str, range(-span_rad, span_rad + 1, res))))
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',', header=head)

def _peak_profile_helper(wlist_all, resolution, fileout):
    """ Helper function for peak_profile_bp_resolution() or peak_profile_wide(). Writes all peak
        profiles to a wiggle file.

    :param wlist_all:
    :param resolution:
    :param fileout:
    """
    with open(fileout + "_bpeaks.wig", 'w') as wigout:
        wlist_all = np.asarray(wlist_all)
        wlist_all = wlist_all[np.lexsort((wlist_all[:, 1].astype(int), wlist_all[:, 0])), :]
        chr_prev = None
        for i in range(wlist_all.shape[0]):
            row = wlist_all[i, :]
            chr_i = row[0]
            if chr_prev != chr_i:
                chr_prev = chr_i
                wigout.write("variableStep\tchrom=%s\n" % chr_i)
            sta_i = int(row[1])
            wlist = row[2:].astype(float)
            for j, x in enumerate(wlist):
                wigout.write("%i\t%0.5f\n" % (sta_i + j * resolution, x))

def per_site_enrchment_ontgt(generator, genome, bamfilein, bamctrlin, fileout, span_rad):
    """
    Function to compute enrichment of a particular factor from ChIP data on a per-site basis.
    :param generator: generator output of msa.target_gen , of the form rs, cut, sen, pam, gui, mis, guide
    :param genome: list of the form [genome_name, path_to_genome.fa]
    :param bamfilein: input bamfile
    :param bamctrlin: control bamfile to use as reference for computing enrichment
    :param fileout: output csv file with cut site
    :param span_rad:
    :return:
    """
    hgsize = get_genome_dict(genome[0])
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    bamctrl = pysam.AlignmentFile(bamctrlin, 'rb')
    wlist, csv_peaks = None, []
    norm_in = bamin.mapped / 1E6
    norm_ctrl = bamctrl.mapped / 1E6
    for rs, cut, sen, pam, gui, mis, guide in generator:
        chr_i = re.split('[:-]', rs)[0]
        sta_i = max(1, cut - span_rad)
        end_i = min(cut + span_rad, hgsize[chr_i])
        rs_i = "%s:%i-%i" % (chr_i, sta_i, end_i)
        cnt_smpl = bamin.count(region=rs_i) / norm_in
        cnt_ctrl = bamctrl.count(region=rs_i) / norm_ctrl
        enrhcmnt = cnt_smpl - cnt_ctrl
        csv_peaks.append([chr_i, cut] + [cnt_smpl, cnt_ctrl, enrhcmnt])
    bamin.close()
    bamctrl.close()
    head = ",".join(["chr", "cut", "Counts EP (RPM)", "Counts Ctrl (RPM)", "Difference EP-NoEP (RPM)"])
    np.savetxt(fileout + "_per-site.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',', header=head)

def per_site_enrchment(list_tgts, genome, bamfilein, bamctrlin, span_rad, fileout):
    """
    Function to compute enrichment of a particular factor from ChIP data on a per-site basis.
    Same as above, but using macs peaks as input.
    :param list_tgts: list of targets
    :param genome: list of the form [genome_name, path_to_genome.fa]
    :param bamfilein: input bamfile
    :param bamctrlin: control bamfile to use as reference for computing enrichment
    :param fileout: output csv file with cut site
    :param span_rad:
    :return:
    """
    hgsize = get_genome_dict(genome[0])
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    bamctrl = pysam.AlignmentFile(bamctrlin, 'rb')
    csv_peaks = []
    norm_in = bamin.mapped / 1E6
    norm_ctrl = bamctrl.mapped / 1E6
    for chr_i, aver_i, span_sta, span_end, fenr_i in list_tgts:
        sta_i = aver_i - span_rad
        end_i = aver_i + span_rad
        if sta_i > 0 and end_i <= hgsize[chr_i]:
            rs_i = "%s:%i-%i" % (chr_i, sta_i, end_i)
            cnt_smpl = bamin.count(region=rs_i) / norm_in
            cnt_ctrl = bamctrl.count(region=rs_i) / norm_ctrl
            d_cnts = cnt_smpl - cnt_ctrl
            norm_d_cnts = d_cnts / float(fenr_i)
            csv_peaks.append([chr_i, aver_i] + [cnt_smpl, cnt_ctrl, d_cnts, norm_d_cnts])
    bamin.close()
    bamctrl.close()
    head = ",".join(["chr", "cut", "Counts EP (RPM)", "Counts Ctrl (RPM)", "Difference EP-NoEP (RPM), "
                                                                           "Normalilzed Difference EP-NoEP (RPM)"])
    np.savetxt(fileout + "_per-site.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',', header=head)

def read_pair_generator(bamfile, region_string=None):
    """
    Generator function that accompanies split_reads_sschip
    :param bamfile: input bam file to split reads
    :param region_string:
    :return:
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bamfile.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def split_reads_sschip(inpbam, outbam):
    """
    Split reads into forward and reverse, for strand-specific ChIP-Seq
    :param inpbam: input bam (RAD51 ChIP-Seq)
    :param outbam: output bam (without extension)
    :return: None
    """
    # Read the input file, define python objects to generate output files
    bamin = pysam.AlignmentFile(inpbam, "rb")
    forbam = pysam.AlignmentFile(outbam + "_for.bam", "wb", template=bamin)
    revbam = pysam.AlignmentFile(outbam + "_rev.bam", "wb", template=bamin)
    for read1, read2 in read_pair_generator(bamin):
        if read1.is_forward:
            forbam.write(read1)
            forbam.write(read2)
        else:
            revbam.write(read1)
            revbam.write(read2)
    bamin.close()
    forbam.close()
    revbam.close()
    # Finally, sort and index.
    pysam.sort("-o", outbam + "_sorted_for.bam", outbam + "_for.bam")
    pysam.sort("-o",  outbam + "_sorted_rev.bam", outbam + "_rev.bam")
    pysam.index(outbam + "_sorted_for.bam")
    pysam.index(outbam + "_sorted_rev.bam")

def get_span_width(generator, genome, f_test, f_ctrl, outpath, w_rad=10000, skip=5000, false_ct=10):
    """ Determine width of 53BP1 or gH2AX enrichment by comparing test sample to negative control
        sample. Extending from the cut site at fixed intervals, enrichment width on either end of
        the cut site is defined to be where there are under 'false_ct' evaluations of negative
        control sample enrichment that is higher than test sample enrichment.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param f_test: test sample BAM file
    :param f_ctrl: negative control BAM file
    :param outpath: path to output BED file (.bed extension omitted)
    :param w_rad: radius of window of enrichment evaluation at each site
    :param skip: number of bases to skip per evaluation of enrichment over control
    :param false_ct: maximum number of times control sample has higher enrichment than test sample
                     for a region to be included in enrichment span width centered at the cut site
    """
    hgsize = get_genome_dict(genome[0])
    outbed = open(outpath + ".bed", 'w')
    outnpy = []
    bam_test, bam_ctrl = pysam.AlignmentFile(f_test, 'rb'), pysam.AlignmentFile(f_ctrl, 'rb')
    cter = 0
    for rs, cut, sen, pam, gui, mis, guide in generator:
        cter += 1
        if cter % 10 == 0:
            print("get_span_width(): Read %i lines of peak file." % cter)
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        index_neg, count_neg, width_neg = 0, 0, 0
        while True:
            index_neg -= skip
            ind_lt_neg, ind_rt_neg = cut + index_neg - w_rad, cut + index_neg + w_rad
            if ind_lt_neg >= 0:
                rs_neg = chr_i + ":" + str(ind_lt_neg) + "-" + str(ind_rt_neg)
                rpm_neg_test = bam_test.count(region=rs_neg) / bam_test.mapped * 1E6
                rpm_neg_ctrl = bam_ctrl.count(region=rs_neg) / bam_ctrl.mapped * 1E6
                if rpm_neg_test <= rpm_neg_ctrl:
                    count_neg += 1
                if count_neg >= false_ct:
                    break
            else:
                break
        index_pos, count_pos, width_pos = 0, 0, 0
        while True:
            index_pos += skip
            ind_lt_pos, ind_rt_pos = cut + index_pos - w_rad, cut + index_pos + w_rad
            if ind_rt_pos <= hgsize[chr_i]:
                rs_pos = chr_i + ":" + str(ind_lt_pos) + "-" + str(ind_rt_pos)
                rpm_pos_test = bam_test.count(region=rs_pos) / bam_test.mapped * 1E6
                rpm_pos_ctrl = bam_ctrl.count(region=rs_pos) / bam_ctrl.mapped * 1E6
                if rpm_pos_test <= rpm_pos_ctrl:
                    count_pos += 1
                if count_pos >= false_ct:
                    break
            else:
                break
        span_rs = chr_i + ":" + str(cut + index_neg) + "-" + str(cut + index_pos)
        enrich_test = bam_test.count(region=span_rs) / bam_test.mapped * 1E6
        enrich_ctrl = bam_ctrl.count(region=span_rs) / bam_ctrl.mapped * 1E6
        bed_1, bed_2, bed_3 = chr_i, str(cut + index_neg), str(cut + index_pos)
        bed_4, bed_5, bed_6 = chr_i + ":" + str(cut), "%0.6f" % (enrich_test - enrich_ctrl), "+"
        bed_7, bed_8 = str(sta_i), str(end_i)
        outbed.write("\t".join([bed_1, bed_2, bed_3, bed_4, bed_5, bed_6, bed_7, bed_8]) + "\n")
        outnpy.append([rs, chr_i, str(cut), bed_2, bed_3, str(int(bed_3) - int(bed_2))])
    header = "region_string, chromosome, cut coordinate, " \
             "start coordinate, end coordinate, span width"
    np.savetxt(outpath + ".csv", np.asarray(outnpy), fmt='%s', delimiter=',', header=header)
    bam_test.close()
    bam_ctrl.close()
    outbed.close()