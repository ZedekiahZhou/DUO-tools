#!/bin/bash
# change log:
#   v1: 2023/11, an integrated pipeline for data cleanning and glori
#   v2: 2023/12/16, supporting Target-seq (--tools, --species, --combine); supportting duplicates not removing (--no_rmdup); 
# shellcheck disable=SC2016
# shellcheck disable=SC2086

ARGS=$(getopt -o 'hi:d:g:p:t:m:' \
  -l 'help,inputfile:,downsample:,species:,thread:,tools:,umi5_len:,umi3_len:,min_len:,fastqc,prx:,untreat,combine,clean_dir:,glori_dir:,only_clean,no_rmdup,gzip,only_glori,m6Am,pipeline:,test' \
  -- "$@")
if [ $? != 0 ] ; then echo "Parse error! Terminating..." >&2 ; exit 1 ; fi
eval set -- "$ARGS"

USAGE='
usage: glori [options]

-h --help		print this help text and exit
-i --inputfile		raw fastq file for cleaning and whole pipeline; clean fastq file for only_glori
-d --downsample N	downsample to N reads
-g --species		genome used (Now hg38 and mm10 supported) or config file indicating information for reference files
-m --min_len N		discard reads that bacame shorter than length N when trimming
-p --thread		threads used
-t --tools		default is "STAR" (which means mapping to genome using STAR and mapping to transcriptome using bowtie), change to bowtie for Target-seq
--umi5_len N		length of bases to removed from 5 prime of the reads
--umi3_len N		length of based to removed from 3 prime of the reads
--fastqc		run fastqc on raw, trimmed and clean fastq file
--prx			prefix of output file name, default is the prefix of input file name
--untreat		run glori in untreat mode (only mapping to sense genome and transcriptome, do not call m6A)
--combine		whether mapping to transcriptome
--clean_dir		output directory for cleaned fastq files, default is ../02_Clean/
--glori_dir		output directory for glori outputs, default is ../03_Sites/
--only_clean		only run data cleanning (trim adpater, deduplication, trim umi)
--no_rmdup		do not remove duplicates
--gzip			compress the output      
--only_glori		only run GLORI (${clean_dir}/${prx}_clean.fq must exits!!)
--m6Am			run GLORI with no filtering in Ratio and P, "-r 0 -p 1.1 -adp 1.1"
--pipeline		run part of GLORI 
--test			only output command to be run without actually running it
'


# parameters -----
while true; do
  case $1 in
    -h|--help) printf "%s\\n" "$USAGE"; exit 2;;
    -i|--inputfile) inputfile=$2; shift 2;;
    -d|--downsample) ds=$2; shift 2;;
    -g|--species) species=$2; shift 2;;
    -m|--min_len) min_len=$2; shift 2;;
    -p|--thread) thread=$2; shift 2;;
    -t|--tools) tools=$2; shift 2;;
    --umi5_len) umi5=$2; shift 2;;
    --umi3_len) umi3=$2; shift 2;;
    --fastqc) fastqc=true; shift;;
    --prx) prx=$2; shift 2;;
    --untreat) untreat=true; shift;;
    --combine) combine=true; shift;;
    --clean_dir) clean_dir=$2; shift 2;;
    --glori_dir) glori_dir=$2; shift 2;;
    --only_clean) no_glori=true; shift;;
    --no_rmdup) no_rmdup=true; shift;;
    --gzip) gzip=true; shift;;
    --only_glori) no_clean=true; shift;;
    --m6Am) m6Am=true; shift;;
    --pipeline) pipeline=$2; shift 2;;
    --test) test=true; shift;;
    --) shift; break ;;
    *) echo "Unrecognized options $1!"; exit 1;;
  esac
done

# 0. set default value for variables
clean_dir=${clean_dir:-'../02_Clean'}
glori_dir=${glori_dir:-'../03_Mapped/'}
mkdir -p ${clean_dir} ${glori_dir}
tools=${tools:-'STAR'}
pipeline=${pipeline:-'All'}
thread=${thread:-'20'}
min_len=${min_len:-'25'}

if [ ! $prx ]; then
  prx=$(basename ${inputfile} | sed -E 's/(_L[1-9]_[12])?.f(ast)?q(.gz)?//')
fi

case $species in
  hg38|mm10) source "/lustre2/chengqiyi_pkuhpc/zhouzhe/toolkit/config/${species}.conf";;
  "") echo "Please specifiy genome used (hg38 or mm10)!";;
  *) source $species;;
esac

# run commands ------
echo $prx

if [ ! $no_clean ]; then
  ## fastqc on raw file
  if [ $fastqc ]; then
    cmd="mkdir -p ../fastqc/raw ../fastqc/clean ../fastqc/trimmed;
    fastqc ${inputfile} --thread ${thread} -q -o ../fastqc/raw"
    echo $cmd | tr ';' '\n' | sed 's/^ //'
    if [ ! ${test} ]; then eval $cmd; fi
  fi

  ## downsample & trim adapter
  if [ $ds ]; then
    cmd="zcat ${inputfile} | seqkit head -n $ds -o ${clean_dir}/${prx}_ds.fq.gz;
    trim_galore --nextseq 20 -j 7 --stringency 1 -e 0.3 --length $((umi5+umi3+min_len)) -o ${clean_dir} ${clean_dir}/${prx}_ds.fq.gz;
    mv ${clean_dir}/${prx}_ds_trimmed.fq.gz ${clean_dir}/${prx}_trimmed.fq.gz"
  else
    tmpfile=$(basename ${inputfile} | sed -E 's/.f(ast)?q(.gz)?/_trimmed.fq.gz/')
    cmd="trim_galore --nextseq 20 -j 7 --stringency 1 -e 0.3 --length $((umi5+umi3+min_len)) -o ${clean_dir} ${inputfile};
    mv ${clean_dir}/${tmpfile} ${clean_dir}/${prx}_trimmed.fq.gz"
  fi
  echo $cmd | tr ';' '\n' | sed 's/^ //'
  if [ ! ${test} ]; then eval $cmd; fi;

  ## fastqc on trimmed file
  if [ $fastqc ]; then
    cmd="fastqc ${clean_dir}/${prx}_trimmed.fq.gz --thread ${thread} -q -o ../fastqc/trimmed"
    echo $cmd | tr ';' '\n' | sed 's/^ //'
    if [ ! ${test} ]; then eval $cmd; fi;
  fi

  ## remove duplicates
  if [ ! ${no_rmdup} ]; then
    if [ ${gzip} ]; then
      cmd="seqkit rmdup -j 10 -s ${clean_dir}/${prx}_trimmed.fq.gz | fastx_trimmer -Q 33 -f $((umi5+1)) -z -o ${clean_dir}/${prx}_clean.fq.gz";
    else
      cmd="seqkit rmdup -j 10 -s ${clean_dir}/${prx}_trimmed.fq.gz | fastx_trimmer -Q 33 -f $((umi5+1)) -o ${clean_dir}/${prx}_clean.fq";
    fi
  else
    if [ ${gzip} ]; then
      cmd="zcat ${clean_dir}/${prx}_trimmed.fq.gz | fastx_trimmer -Q 33 -f $((umi5+1)) -z -o ${clean_dir}/${prx}_clean.fq.gz";
    else
      cmd="zcat ${clean_dir}/${prx}_trimmed.fq.gz | fastx_trimmer -Q 33 -f $((umi5+1)) -o ${clean_dir}/${prx}_clean.fq";
    fi
  fi
  echo $cmd
  if [ ! ${test} ]; then eval $cmd; fi;

  ## fastqc on clean file
  if [ $fastqc ]; then
    cmd="fastqc ${clean_dir}/${prx}_clean.fq* --thread ${thread} -q -o ../fastqc/clean"
    echo $cmd | tr ';' '\n' | sed 's/^ //'
    if [ ! ${test} ]; then eval $cmd; fi;
  fi
fi

if [ ! $no_glori ]; then
  if [ ! $no_clean ]; then 
    clean_file=${clean_dir}/${prx}_clean.fq
  else
    clean_file=${inputfile}
  fi

  cmd_prx="python ${tooldir}/run_GLORI.py -i $tooldir -q ${clean_file} -T $thread -f ${genome} -t ${tools} --pipeline ${pipeline} -pre ${prx} -o ${glori_dir}"
  if [ ${untreat} ]; then
    if [ ${combine} ]; then
      cmd="${cmd_prx} -Tf ${TfGenome} -a $anno --combine --untreated"
    else
      cmd="${cmd_prx} --untreated"
    fi
  else
    if [ ${m6Am} ]; then 
      add_pra='-c 1 -C 0 -r 0 -p 1.1 -adp 1.1 -s 0'
    else
      add_pra=''
    fi

    if [ ${combine} ]; then
      cmd="${cmd_prx} -f2 ${genome2} -rvs ${rvsgenome} -Tf ${TfGenome} -a $anno -b $baseanno --combine --rvs_fac ${add_pra}"
    else
      cmd="${cmd_prx} -f2 ${genome2} ${add_pra}"
    fi
  fi
  echo $cmd
  if [ ! ${test} ]; then eval $cmd; fi;
fi



