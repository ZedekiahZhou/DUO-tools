"""
Author: Zhe Zhou, Peking University, Yi lab
Date: May 24, 2024
Email: zzhou24@pku.edu.cn
Program: This program  is used for DUO-seq analysis
Version: 1.0.0
"""

import argparse, subprocess, re, os, time, sys
from time import strftime

parser = argparse.ArgumentParser(description="DUO-tools for detecting m6Am sites from DUO-seq data")

group_global = parser.add_argument_group("Global")
group_global.add_argument("--mode", dest="mode", default="all", choices=["all", "preprocessing", "mapping", "call_TSS", "call_m6Am"], 
                          help='Run which part of the pipeline: all, preprocessing, mapping, call_TSS or call_m6Am; separated by comma if multiple part supplied')
group_global.add_argument("--raw_fq", dest="raw_fq", help="raw fastq(.gz) file")
group_global.add_argument("--clean_fq", dest="clean_fq", help="clean fastq(.gz) file for mapping")
group_global.add_argument("--bam", dest="bam", help="merged sorted bam files")
group_global.add_argument("--prx", dest="prx", help="output file prefix")
group_global.add_argument("--DUOdir", dest="DUOdir", default=os.path.dirname(__file__)+"/", help="directory of DUO-tools")
group_global.add_argument("-o", "--outdir", dest="outdir", default="./res/", help="output directory, default is ./")
group_global.add_argument("-p", "--threads", dest="threads", default=20, help="threads used, default is 20")

group_preprocessing = parser.add_argument_group("Preprocessing")
group_preprocessing.add_argument("--umi5", dest="umi5", default=0, help="length of bases to removed from 5 prime of the reads, default is 0")
group_preprocessing.add_argument("--umi3", dest="umi3", default=0, help="length of bases to removed from 5 prime of the reads, default is 0")
group_preprocessing.add_argument("-m", "--min_len", dest="min_len", default=25, help="discard reads that bacame shorter than min_len before mapping")
group_preprocessing.add_argument("--fastqc", dest="fastqc", default=True, help="whether to run fastqc, default is True")
group_preprocessing.add_argument("--tag_seq", dest="tag_seq", default="TGACGCTGCCGACGATC", 
                                 help="tag sequence ligated to 5' ends of TSS, default is TGACGCTGCCGACGATC")

group_mapping = parser.add_argument_group("Mapping")
group_mapping.add_argument("-f", "--reference", nargs="?", help="Index file for the plus strand of the genome")
group_mapping.add_argument("-f2", "--reference2", nargs="?", help="Index file for the unchanged genome")
group_mapping.add_argument("-rvs", "--rvsref", nargs="?", help = "Index file for the minus strand of the genome")
group_mapping.add_argument("-Tf", "--transref", nargs="?", help="Index file for the minus strand of the transcriptome")
group_mapping.add_argument("-a", "--anno", nargs="?", help="Annotation file within exons")
    
group_m6Am = parser.add_argument_group("Call m6Am")
group_m6Am.add_argument("-ta", "--tssanno", help="Annotation of TSS range")

args = parser.parse_args()

def fun_pre(raw_fq, prx, args):
    """"""
    print(prx, flush=True)
    print("\n[%s] Preprocessing ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)

    if args.fastqc:
        #subprocess.call("mkdir -p " + args.outdir +"/fastqc/raw " + args.outdir + "/fastqc/clean", shell=True)
        cmd="fastqc "+ raw_fq + " --thread " + str(args.threads) + " -q -o " + args.outdir +"/fastqc/raw"
        print(cmd, flush=True)
        #subprocess.call(cmd, shell=True)

    clean_dir=args.outdir + "/02_Clean/"

    # trim adapter
    cutoff_len1=args.umi5+args.umi3
    cmd="trim_galore -q 20 -j 7 --stringency 1 -e 0.3 --length " + str(cutoff_len1) + " -o " + clean_dir + " " + raw_fq
    print(cmd, flush=True)
    #subprocess.call(cmd, shell=True)

    tmpfile=clean_dir + re.match(".+/([^/]+).f(ast)?q(.gz)?$", raw_fq).group(1) + "_trimmed.fq.gz"
    trimmed_fq=clean_dir + prx + "_trimmed.fq.gz"
    cmd="mv " + tmpfile + " " + trimmed_fq
    print(cmd, flush=True)
    #subprocess.call(cmd, shell=True)

    # remove duplicates and umi
    rmdup_fq=clean_dir + prx + "_rmdup.fq.gz"
    cmd="seqkit rmdup -j 10 -s " + trimmed_fq + " | fastx_trimmer -Q 33 -f " + str(args.umi5+1) + " -z -o " + rmdup_fq
    print(cmd, flush=True)
    #subprocess.call(cmd, shell=True)

    # remove tag
    clean_fq=clean_dir + prx + "_clean.fq"
    tag_info=clean_dir + prx + "_rmtag.info"
    cmd='cutadapt -j 0 -g "' + args.tag_seq + ';rightmost" -m ' + str(args.min_len) + ' -O ' + str(len(args.tag_seq)) + \
        ' -e 0.2 -o ' + clean_fq + ' --info-file ' + tag_info + ' --discard-untrimmed ' + rmdup_fq
    print(cmd, flush=True)
    subprocess.call(cmd, shell=True)
    
    if args.fastqc:
        cmd="fastqc "+ clean_fq + " --thread " + str(args.threads) + " -q -o " + args.outdir +"/fastqc/clean"
        print(cmd, flush=True)
        subprocess.call(cmd, shell=True) 


def fun_mapping(clean_fq, prx, args):
    site_dir=args.outdir + "/03_Sites/"
    cmd_prx="python " + args.DUOdir + "/Mapping/Run_GLORI_Mapping.py -q " + clean_fq + " -T " + str(args.threads) + \
        " -f " + args.reference + " -pre " + prx + " -o " + site_dir
    cmd=cmd_prx + " -f2 " + args.reference2 + " -rvs " + args.rvsref+ " -Tf " + args.transref + " -a " + args.anno + " --combine --rvs_fac"
    print(cmd, flush=True)
    subprocess.call(cmd, shell=True)


def fun_m6Am(bam, prx, args):
    site_dir=args.outdir + "/03_Sites/"
    
    # call TSS
    ftss=site_dir + prx + "_TSS_raw.bed"
    cmd="python " + args.DUOdir + "/Pileup/get_TSS.py -r " + args.reference2 + " -b " + bam + " -o " + ftss
    print(cmd, flush=True)
    subprocess.call(cmd, shell=True)

    # annotate TSS
    ftss_anno=site_dir + prx + "_TSS_raw.bed.annotated"
    cmd="bedtools intersect -a " + ftss + " -b " + args.tssanno + " -s -wa -wb > " + ftss_anno
    print(cmd, flush=True)
    subprocess.call(cmd, shell=True)

    cmd="python " + args.DUOdir + "/Pileup/anno_TSS.py " + ftss_anno
    print(cmd, flush=True)
    subprocess.call(cmd, shell=True)

    # 筛选TSS
    
    # pile ATCG
    ftss_clean=site_dir + prx + "_TSS.bed"
    fAGcount=site_dir + prx + "_AGcount.tsv"
    cmd="python " + args.DUOdir + "/Pileup/pileup_reads5p.py -r " + args.reference2 + " -l " + ftss_clean + " -o " + fAGcount + " -b " + bam

    # call m6Am



def main():
    global args
    mode = args.mode.split(',')
    os.makedirs(args.outdir+"/02_Clean/", exist_ok=True)
    os.makedirs(args.outdir+"/03_Sites/", exist_ok=True)
    os.makedirs(args.outdir+"/intermediate/", exist_ok=True)

    if "all" in mode:
        mode = ["preprocessing", "mapping", "call_m6Am"]
    prx=None

    if "preprocessing" in mode:
        if args.raw_fq is None:
            raise ValueError("Raw fastq files must be provided if run preprocessing!")
        if args.prx is None:
            prx=re.match("(.*/)?([^/]+)_L[1-9]+_(R)?[12].f(ast)?q(.gz)?$", args.raw_fq).group(2)
        else:
            prx=args.prx
        fun_pre(args.raw_fq, prx, args)

    if "mapping" in mode:
        if prx is None:  # begin with mapping step
            if args.clean_fq is None:
                raise ValueError("Clean fastq files must be provided if beginning with the mapping step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_clean.fq$", args.clean_fq).group(2)
                else:
                    prx=args.prx
        else:
            clean_fq=args.outdir + "/02_Clean/" + prx + "_clean.fq"
            fun_mapping(clean_fq, prx, args)

    if "call_m6Am" in mode:
        if prx is None:  # begin with m6Am calling step
            if args.bam is None:
                raise ValueError("Bam files must be provided if beginning with the m6Am calling step!")
            else:
                if args.prx is None:
                    prx=re.match("(.*/)?([^/]+)_merged.sorted.bam$", args.bam).group(2)
                else:
                    prx=args.prx
        else:
            bam=args.outdir + "/03_Sites/" + prx + "_merged.sorted.bam"
            fun_m6Am(bam, prx, args)


if __name__ == '__main__':
    main()
