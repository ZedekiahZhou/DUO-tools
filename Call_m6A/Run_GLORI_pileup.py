"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to run GLORI-tools"""
"""Input: [.fastq]"""

# version information
__version__ = "0.1.2"

import os, sys
import argparse
import time
import subprocess
import os
from time import strftime


def run_command(finalbam, Threads):
    f_pileup = outputprefix + ".pileup"  # pileup file

    print("\n[%s] Pileup ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("python "+DUOdir+"pileup_genome_multiprocessing.py -P " + Threads + " -f " + genome + " -i " + finalbam + " -o " + f_pileup)
    subprocess.call("python "+DUOdir+"pileup_genome_multiprocessing.py -P "+Threads+" -f "+genome+" -i "+ finalbam +" -o "+f_pileup,shell=True)
    print("python "+DUOdir+"get_referbase.py -input " + f_pileup + " -referFa "+ genome2 + " -outname_prx "\
                    + outputdir +'/'+prx)
    subprocess.call("python "+DUOdir+"get_referbase.py -input " + f_pileup + " -referFa "+ genome2 + " -outname_prx "\
                    + outputdir +'/'+prx, shell=True)
    subprocess.call("rm -f " + f_pileup, shell=True)

    print("\n[%s] Done! ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GLORI-tools for pileup")
    parser.add_argument("-T", "--Threads", nargs="?", type=str, default='1',help="Used threads")

    group_required = parser.add_argument_group("Required")
    group_required.add_argument("--bam", nargs="?", type=str, default=sys.stdin, help="Merged sorted bam file from Run_GLORI_Mapping.py")
    group_required.add_argument("-f", "--reference", nargs="?", type=str, default=sys.stdin, help="Index file for the plus strand of the genome")
    group_required.add_argument("-f2", "--reference2", nargs="?", type=str, default=sys.stdin, help="Index file for the unchanged genome")

    group_output = parser.add_argument_group("Output (optional)")
    group_output.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default', help="Output file prefix")
    group_output.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="Output directory")
    

    options = parser.parse_args()
    global genome,genome2,outputdir,Threads,prx,outputprefix
    DUOdir = os.path.dirname(__file__)+"/"
    genome = options.reference
    genome2 = options.reference2
    outputdir = options.outputdir
    Threads = options.Threads
    prx = options.outname_prefix
    
    
    outputprefix = outputdir + "/" + prx
    run_command(options.bam, Threads)


