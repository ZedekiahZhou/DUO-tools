#!/bin/bash

source /lustre2/chengqiyi_pkuhpc/zhouzhe/software/DUO/config/mm10.conf
python ${DUOdir}/DUO.py --raw_fq ${file} -o ../ -f ${genome} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a ${anno} -ta ${tssanno}