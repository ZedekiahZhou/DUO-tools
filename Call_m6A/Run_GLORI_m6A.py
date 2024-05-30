"""Cong Liu Yi_lab, Peking University"""
"""Feb, 2021"""
"""Email: liucong-1112@pku.edu.cn"""
"""Usage: This program is used to run GLORI-tools"""
"""Input: [.fastq]"""

# version information
__version__ = "0.1.2"

"""
History:
    mod_v1.2: 
        1) add parameters to control pipeline
    mod_v1.1: 240326
        1) use awk to do A-G change
        2) make '-g' parameter usable for m6A_caller.py
    mod_v1.0: 
        1) Add comments; 
        2) not delete Important tmp files!!
        3) When coverge is < 15 for specifi segment (eg. chrM), pandas throw an error when output *.totalCR.txt, add an if-else to advoid.
        4) reformat log for timming
"""

import pandas as pd
import os, sys
import argparse
import time
import subprocess
import os
import multiprocessing
from time import strftime
import re

def get_sites(total_mpi, chr):  
    """format pileup file and call m6A sites for each chromosome"""
    print(chr)
    # time.sleep(1)
    chr2 = chr.split("_AG_converted")[0]
    baseanno_chr = baseanno + "." + chr2
    file_mpi = total_mpi + "." + chr2
    file_format = outputprefix + ".referbase.mpi.formatted.txt" + "." + chr2
    file_CR = outputprefix + ".CR.txt" + "." + chr2
    file_sites = outputprefix + ".callsites" + "." + chr2
    print("awk \'$1==\"\'\"" + chr + "\"\'\"\' " + total_mpi + " > " + file_mpi)
    subprocess.call("awk \'$1==\"\'\"" + chr + "\"\'\"\' " + total_mpi + " > " + file_mpi,shell = True)
    if baseanno != 'None':
        if not os.path.exists(baseanno_chr):
            print("awk \'$1==\"\'\"" + chr2 + "\"\'\"\' " + baseanno + " > " + baseanno_chr)
            subprocess.call("awk \'$1==\"\'\"" + chr2 + "\"\'\"\' " + baseanno + " > " + baseanno_chr,shell = True)
        print("python "+DUOdir+"m6A_pileup_formatter.py --db "+ baseanno_chr + " -i " + file_mpi + " -o "+ file_format +" --CR "+ file_CR)
        subprocess.call("python "+DUOdir+"m6A_pileup_formatter.py --db "+ baseanno_chr + " -i " + file_mpi + " -o "+file_format +" --CR "+ file_CR,shell=True)
    else:
        print("python " + DUOdir + "m6A_pileup_formatter.py " + " -i " + file_mpi + " -o " + file_format + " --CR " + file_CR)
        subprocess.call("python " + DUOdir + "m6A_pileup_formatter.py " + " -i " + file_mpi + " -o " + file_format + " --CR " + file_CR,
            shell=True)
    print("python "+DUOdir+"m6A_caller.py -i "+file_format + " -o " + file_sites +" -c " + Cov +" -C "+ Counts + " -r "\
                    + minRatio +" -p "+pvalue+" -s "+multiAratio+" -R "+AGRatio +" --cutoff " + Acutoffs +" --CR "+ \
                    background+" --method " + statmethod + " -g " + geneCR)
    subprocess.call("python "+DUOdir+"m6A_caller.py -i "+file_format + " -o " + file_sites +" -c " + Cov +" -C "+ Counts + " -r "\
                    + minRatio +" -p "+pvalue+" -s "+multiAratio+" -R "+AGRatio +" --cutoff " + Acutoffs +" --CR "+ \
                    background+" --method " + statmethod + " -g " + geneCR,shell=True)


def run_command(total_mpi, Threads):
    chr_file = outputprefix + "_chrlist"
    final_sites1 = outputprefix + ".totalm6A.txt"  # sites with raw P < pcutoff
    final_format = outputprefix + ".totalformat.txt"  # pooled formatted pileup file
    final_sites2 = outputprefix + ".totalm6A.FDR"  # FDR filtered m6A sites (actually nearly no sites were removed duo to wrong FCR calculation!!)

    # call sites
    print("\n[%s] Call sites ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("cut -f 1 "+total_mpi + " | sort -u > " + chr_file)
    subprocess.call("cut -f 1 "+total_mpi + " | sort -u > " + chr_file, shell=True)
    chr_list = sorted([r1.strip().split("\t")[0] for r1 in open(chr_file).readlines()])
    #
    fac_chr = True
    index = 0

    # format pileup file and call m6A sites for each chromosome
    while fac_chr:
        chr_prx = chr_list[0][:3]
        detected_chr = [i.split('.')[-1]+"_AG_converted" for i in os.listdir(outputdir) if re.search(prx + ".referbase.mpi.formatted.txt." + chr_prx,i) and not re.search('tmp',i)]
        print("*************Detected",detected_chr,len(detected_chr))
        if len(detected_chr) == len(chr_list):
            print(index)
            fac_chr = False
            break
        print(index,'lalalalala')
        if index >= 12:
            print("Erro!!! break multiprocessing")
            fac_chr = False
        else:
            index += 1
            differ_chr = [i for i in chr_list if i not in detected_chr]
            print("**************No_Detected",differ_chr)
            multiprocessing.freeze_support()
            pool = multiprocessing.Pool(int(Threads))
            try:
                for chr in differ_chr:
                    pool.apply_async(func=get_sites, args=(total_mpi, chr,))
                pool.close()
                pool.join()
            finally:
                pool.terminate()

    # Combine pileup file
    print("\n[%s] Conbine results ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    # subprocess.call("rm -f " + total_mpi, shell=True)  # do not remove *.referbase.mpi
    subprocess.call("mkdir -p "+outputdir+"/intermediate", shell=True)
    subprocess.call("mv "+ total_mpi + " " + outputdir + "/intermediate", shell=True)
    
    print("cat " + outputprefix + ".referbase.mpi.formatted.txt" + ".* | sed '/#/d' > " + final_format)
    subprocess.call("cat " + outputprefix + ".referbase.mpi.formatted.txt" + ".* | sed '/#/d' > " + final_format,
                    shell=True)
    
    # Combing A-to-G conversion files
    print("**************Combing A-to-G conversion files******************")
    pd_CR = pd.DataFrame()
    final_CR = outputprefix + ".totalCR.txt"

    for CR_c in chr_list:
        chr_x = CR_c.split("_AG_converted")[0]
        file = outputprefix+".CR.txt."+chr_x
        CR1 = pd.read_csv(file, sep="\t", names=['SA', 'Totalcovered_reads', 'Remained A reads', 'Non-A-to-G ratio','Mapped_area'])

        if CR1.shape[0] > 0:  # CR.txt files for some segment may be NULL!!
            CR2 = CR1[~CR1['SA'].isin(['#90%', '#75%', '#50%', '#25%', '#10%', '#Median', '#Mean'])]
            CR2.iloc[0, 0] += "_" + chr_x
            pd_CR = pd.concat([pd_CR, CR2])

    # whole transcriptome
    pd_chr = pd_CR[pd_CR.SA.str.contains("#ALL_", case=True)]
    pd_tmp = pd_chr.sum().to_frame().T
    pd_tmp["SA"] = "#ALL_genome"
    pd_tmp["Non-A-to-G ratio"] = pd_tmp["Remained A reads"]/pd_tmp["Totalcovered_reads"]
    pd_CR = pd.concat([pd_tmp, pd_CR])
    pd_CR['A-to-G_ratio'] = 1 - pd_CR['Non-A-to-G ratio']
    pd_CR = pd_CR[['SA', 'A-to-G_ratio', 'Totalcovered_reads', 'Remained A reads', 'Mapped_area']]
    pd_CR.to_csv(final_CR,sep="\t",index=False)


    # combine m6A sites
    print("cat " + outputprefix + ".callsites" + "*." + str(Acutoffs) + ".txt > " + final_sites1)
    subprocess.call("cat " + outputprefix + ".callsites" + "*." + str(Acutoffs) + ".txt > " + final_sites1, shell=True)
    print("rm -f " + outputprefix + "*." + chr_prx + "* ")
    subprocess.call("rm -f " + outputprefix + "*." + chr_prx + "* ", shell=True)
    #subprocess.call("rm -f " + chr_file, shell=True)
    #subprocess.call("mv " + outputprefix + "*." + chr_prx + "* " + outputdir + "/intermediate", shell=True)
    subprocess.call("mv " + chr_file + " " + outputdir + "/intermediate", shell=True)

    # FDR filter
    print("\n[%s] FDR filtering ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    print("python " + DUOdir + "m6A_caller_FDRfilter.py -i " + final_sites1 + ' -o ' + final_sites2\
                    + ' -adp ' + adjP)
    subprocess.call("python " + DUOdir + "m6A_caller_FDRfilter.py -i " + final_sites1 + ' -o ' + final_sites2\
                    + ' -adp ' + adjP, shell=True)
    #subprocess.call("rm -f " + final_sites1, shell=True)
    subprocess.call("mv " + final_sites1 + " " + outputdir + "/intermediate", shell=True)

    print("\n[%s] Done! ========" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GLORI-tools for detecting m6A sites from high-throughput sequencing data")

    group_required = parser.add_argument_group("Required")
    group_required.add_argument("--mpi", nargs="?", type=str, default=sys.stdin, help="Pileup file with reference base")

    group_output = parser.add_argument_group("Output (optional)")
    group_output.add_argument("-pre", "--outname_prefix", nargs="?", type=str, default='default', help="Output file prefix")
    group_output.add_argument("-o", "--outputdir", nargs="?", type=str, default='./', help="Output directory")

    group_mappingfilter = parser.add_argument_group("Mapping conditions (optional)")
    group_mappingfilter.add_argument("-b", "--baseanno", nargs="?", type=str, default='None', help="Annotations at single-base resolution")
    group_mappingfilter.add_argument("-T", "--Threads", nargs="?", type=str, default='1',help="Used threads")

    # Filter
    group_site = parser.add_argument_group("m6A filter (optional)")
    group_site.add_argument("-c", "--coverage", dest="coverage", default='15', type=str, help="A+G coverage")
    group_site.add_argument("-C", "--count", dest="count", default='5', type=str,help="A coverage")
    group_site.add_argument("-r", "--ratio", dest="ratio", default='0.1', type=str, help="m6A level or A rate")
    group_site.add_argument("-p", "--pvalue", dest="pvalue", default='0.005', type=str, help="P-value cutoff for statistical test")
    group_site.add_argument("-adp", "--adjustpvalue", dest="adjustpvalue", default='0.005', type=str, help="Cutoff for FDR-adjusted p-value")
    group_site.add_argument("-s", "--signal", dest="signal", default='0.8', type=str,
                            help="signal ratio, The reads coverage (below the A-cutoff) divided by the total coverage \
                            is used to exclude false positives located in regions resistant to nitrite treatment")
    group_site.add_argument("-R", "--var_ratio", dest="var_ratio", default='0.8', type=str,
                            help="The A+G coverage divided by the total coverage is used to exclude mapping errors")
    group_site.add_argument("-g", "--gene_CR", dest="gene_CR", default='0.2', type=str,
                            help="Any m6A sites within a gene with an A-to-G conversion rate below the cutoff will be discarded")
    group_site.add_argument("-N", "--AG", dest="AG_number", default=0, type=int,
                            help="Any m6A sites within a gene with an A+G coverage below the cutoff will be discarded")


    # Statistics
    group_stat = parser.add_argument_group("Statistic method (optional)")
    group_stat.add_argument("--method", dest="method", default="binomial", choices=['binomial', 'poisson'],
                            help="Statistical method to test the significance of m6A.")
    group_stat.add_argument("--CR", dest="conversion_rate", default="gene", choices=['gene', 'overall'],
                            help="The control conversion rate used to establish statistical models")
    group_stat.add_argument("--NA", dest="non_anno", default="ELSE",
                            choices=['ELSE', 'Median', 'Mean', 'ALL', 'discard'],
                            help="Determining which CR to use for the sites that are not annotated to any genes")
    group_stat.add_argument("--cutoff", dest="A_cutoffs", default="3", choices=["1","2","3","4","5","6","7","8","9","10","15","20","None"],
                            help="A-cutoffs, The cutoff for the minimum number of remained A bases in each read")
    

    options = parser.parse_args()
    global outputdir,Threads,baseanno,prx,outputprefix,Mapping,Pileup,Callsite
    global Cov,Counts,minRatio,pvalue,adjP,multiAratio,AGRatio,geneCR,AG_number_gene,statmethod,background,NAbackground,Acutoffs
    DUOdir = os.path.dirname(__file__)+"/"
    outputdir = options.outputdir
    Threads = options.Threads
    baseanno = options.baseanno
    prx = options.outname_prefix
    Cov = options.coverage
    Counts = options.count
    minRatio = options.ratio
    pvalue = options.pvalue
    adjP=options.adjustpvalue
    multiAratio = options.signal
    AGRatio = options.var_ratio
    geneCR = options.gene_CR
    AG_number_gene = options.AG_number
    statmethod = options.method
    background = options.conversion_rate
    NAbackground = options.non_anno
    Acutoffs = options.A_cutoffs
    

    outputprefix = outputdir + "/" + prx

    if baseanno == 'None':
        background = 'overall'
    
    run_command(options.mpi, Threads)


