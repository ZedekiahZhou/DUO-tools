#!/bin/bash
# shellcheck disable=SC2154

source /lustre2/chengqiyi_pkuhpc/zhouzhe/software/DUO/config/mm10.conf

# run in m6Am mode
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6Am --module preprocessing,mapping --umi5 3 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}
python ${DUOdir}/DUO.py --bam ${bam_file} -o ../ --mode m6Am --module call_m6Am -ds 14000000 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}

# m6Am untreated
python ${DUOdir}/DUO.py --raw_fq ${file} -o ../ --mode m6Am --untreated \
    -f ${genome} -Tf ${TfGenome} -a ${anno} --test


# run in m6A mode
python ${DUOdir}/DUO.py --raw_fq ${fastq_file} -o ../ --mode m6A --module preprocessing,mapping --umi5 3 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno}
python ${DUOdir}/DUO.py --bam ${bam_file} -o ../ --mode m6A --module call_m6A -ds 14000000 \
    -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno}

# m6A untreated
python ${DUOdir}/DUO.py --raw_fq ${file} -o ../ --mode m6A --untreated \
    --umi5 3 -f ${genome} -Tf ${TfGenome} -a ${anno}