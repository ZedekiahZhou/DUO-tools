#!/bin/bash

species="hg38"
genome2=/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/GLORI/hg38/genome/hg38_chr_only.fa

echo -e "\n[ $(date) ]: MetaPlotR --------"
source activate metaPlotR
metaPlotR_dir=/lustre2/chengqiyi_pkuhpc/zhouzhe/software/metaPlotR/
annot_prx=/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/anno_metaPlotR/${species}/${species}

cat $1 | awk -v OFS='\t' 'NR>1 {print $1, $2-1, $2, $4, 0, $3}' > "$1".bed
perl ${metaPlotR_dir}/annotate_bed_file.pl --bed "$1".bed --bed2 ${annot_prx}_annot.sorted.bed > "$1"_annot.bed
perl ${metaPlotR_dir}/rel_and_abs_dist_calc.pl --bed "$1"_annot.bed --regions ${annot_prx}_region_sizes.txt > "$1"_dist.measures.txt
rm -rf "$1".bed "$1"_annot.bed
source activate zzhou_bio2

# 03. get 5 bases around m6A sites
echo -e "\n[ $(date) ]: Get 5 bases around sites --------"
cat $1 | awk -v OFS='\t' 'NR!=1 {print $1, $2-3, $2+2, $4, 0, $3}'  > "$1"_len5.bed
bedtools getfasta -s -fi ${genome2} -bed "$1"_len5.bed -bedOut > "$1"_len5_motif.bed
rm -rf "$1"_len5.bed