#!/bin/bash
# shellcheck disable=SC2154

source /lustre2/chengqiyi_pkuhpc/zhouzhe/software/DUO/config/mm10.conf

# I. m6Am
# single command
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6Am \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}

# downsample to N reads after mapping
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6Am --module preprocessing,mapping \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}
python ${DUOdir}/DUO.py --bam ${bam_file} -o ../ --mode m6Am --module call_m6Am,QC -ds 14000000 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}

# II. m6A
# single command (takara V2)
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6A --umi5 3 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ba ${baseanno}

# downsample to N reads after mapping
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6A --module preprocessing,mapping --umi5 3 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno}
python ${DUOdir}/DUO.py --bam ${bam_file} -o ../ --mode m6A --module call_m6A,QC -ds 14000000 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ba ${baseanno}