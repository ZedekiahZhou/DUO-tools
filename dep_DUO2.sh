#!/bin/bash

ARGS=$(getopt -o 'hr:b:p:' -l 'help,ref:,bam:,thread:,tx_anno:' -- "$@")
if [ $? != 0 ]; then echo "Parse error! Terminating..." >&2; exit 1; fi
eval set -- "$ARGS"

USAGE='
usage: DUO [options]

-h --help           print this help text and exit
-r --ref            reference genome files, seperated by comma if multiple files are provided
-b --bam            use merged sorted bam files as input, seperated by comma if multiple files are provided
-p --thread         number of threads to used
--tx_anno           annotation file of transcripts, from upstream 100nt to transcripts end
'

# parameters ------
while true; do
    case $1 in
        -h|--help) printf "%s\\n" "$USAGE"; exit 2;;
        -r|--ref) refs=$2; shift 2;;
        -b|--bam) fbams=$2; shift 2;;
        -p|--thread) thread=$2; shift 2;;
        --tx_anno) tx_anno=$(echo $2 | tr "," " "); shift 2;;
        -g|--species) species=$2; shift 2;;
        --) shift; break ;;
        *) echo "Unrecognized options $1!"; exit 1;;
    esac
done


DUOdir="/lustre2/chengqiyi_pkuhpc/zhouzhe/software/DUO"
#f_anno="/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/GRCh38/anno/gencode.v43.annotation.TSS.u100dInf.bed"
#refs="/lustre2/chengqiyi_pkuhpc/zhouzhe/genome/GLORI/hg38/genome/hg38_chr_only.fa"
clean_file=""
tools=""
Mapping_dir="../03_Mapping/"
TSS_dir="../04_TSS/"

thread=${thread:-"20"}


case $species in
  hg38|mm10) source "${DUOdir}/config/${species}.conf";;
  "") echo "Please specifiy genome used (hg38 or mm10)!";;
  *) source $species;;
esac



prx=$(basename ${fbams} | sed 's/_merged.+//')
echo $prx

thread=20
ftss=${prx}_TSS.bed

# mapping
python ${DUOdir}/Mapping/Run_GLORI_Mapping.py 
cmd_prx="python ${DUOdir}/Mapping/Run_GLORI_Mapping.py -q ${clean_file} -T ${thread} -f ${AGconverted_genome} -t ${tools} -pre ${prx} -o ${Mapping_dir}"
  if [ ${untreat} ]; then
    if [ ${combine} ]; then
      cmd="${cmd_prx} -Tf ${TfGenome} -a $anno --combine --untreated"
    else
      cmd="${cmd_prx} --untreated"
    fi
  else
    if [ ${combine} ]; then
      cmd="${cmd_prx} -f2 ${original_genome} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno --combine --rvs_fac"
    else
      cmd="${cmd_prx} -f2 ${original_genome} ${add_pra}"
    fi
  fi
  echo $cmd
  if [ ! ${test} ]; then eval $cmd; fi;

# pileup
python ${DUOdir}/Pileup/get_TSS.py -r ${refs} -b ${fbams} -o ${TSS_dir}/${ftss}
bedtools intersect -a ${TSS_dir}/${ftss} -b ${tx_anno} -s -wa -wb > ${TSS_dir}/${ftss}.annotated
python ${DUOdir}/Pileup/anno_TSS.py ${TSS_dir}/${ftss}.annotated

