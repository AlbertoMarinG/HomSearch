#!/bin/bash
# Bash file used to extract genome-wide insulation profiles from juicer .hic files
conda activate cworld
mat2insul="PATH-TO-CWORLD-PERL-FILE (e.g. /home/Codes/cworld-dekker/scripts/perl/matrix2insulation.pl)"
jcrtls="PATH-TO-JUICER-TOOLS (e.g. /home/Codes/juicer/CPU/juicer_tools.jar)"
hicfile="PATH-TO-HIC-FILE (juicer output) (e.g. /home/Data/Hi-C/RAD21_HiC/BR1/neg_no_30.hic)"
chrom_sizes="PATH-TO-CHROM-SIZES-FILE (e.g. /home/Databases/Human_Ref_Gen/GRCh38_noalt_as/hg38.chrom.sizes)"
gen="GENOME-NAME (e.g. hg38)"
binsize="25000"
res=`echo $binsize/1000 | bc`
chr_list=$(egrep -v '(_|Y|M)' $chrom_sizes | cut -f1)
# For each chromosome: 1. Dump matrix from .hic file using juicer, 2. convert to input file needed for cworld, 3. run cworlds' matrix2insulation script
for i in ${chr_list[*]}; do
	outjcr="$i"_"$res"kb.txt
	mat_file=""$i"_"$res"kb_cworld.txt"
	out_ins=""$i"_"$res"kb"
	echo ${i:3}
	java -jar $jcrtls dump -d observed None $hicfile ${i:3} ${i:3} bp $binsize > $outjcr
	mat_file=""$i"_"$res"kb_cworld.txt"
        mat_size=`awk 'NR==1 {print NF}' $outjcr`
	awk '{print NR"|"genome"|"chrom":"init"-"init+binsize-1, $0; init=init+binsize}' OFS="\t" genome=$gen chrom=${i:3} init=1 binsize=$binsize $outjcr > tmp.dat
        cut -f 1 tmp.dat | awk 'BEGIN{ORS="\t"}{print}' | awk '{print bin"x"bin, $0}' OFS="\t" bin=$mat_size > header.dat
        cat header.dat tmp.dat > $mat_file
        rm header.dat tmp.dat
	perl $mat2insul -i $mat_file -o $out_ins
done
# Merge insulation files from all chromosomes into a single file:
outfile=inscore_all.dat
rm $outfile
touch $outfile
for i in ${chr_list[*]}; do
        echo $i
        insfile=""$i"_25kb--is500000--nt0--ids250000--ss0--immean.insulation"
        awk -v OFS='\t' '/^[^#]/ {print chr, $2, $3, $8, $9, $10}' chr=$i $insfile | grep -v "rawInsulationScore" | sed -e 's/NA/0.0/g' > tmp.dat
        cat $outfile tmp.dat > tmp2.dat
        mv tmp2.dat $outfile
        rm tmp.dat
done
~       
