#-------------Dependencies--------------------------

##check if installed and if not "freebayes", "java" and "have snpEff in directory" 
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


import os

import subprocess
import sys
import argparse

try:
    import pandas as pd
except ModuleNotFoundError:
    install("pandas")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from multiprocessing import Pool
import matplotlib.backends.backend_pdf
try:
    import allel
except ModuleNotFoundError:
    install("scikit-allel")
import itertools
import json


def fastq_files(file):

    if os.listdir(file)[0].endswith("fastq.gz"):
        file = file
    elif os.path.isdir(file + "/" + os.listdir(file)[0]): 
        #file = file + "/" + os.listdir(file)[0]
        file = fastq_files(file + "/" + os.listdir(file)[0])
    return file

def variant_calling(out_file, sample , ref):
    print("Performing Variant Calling of Sample {}".format(sample))
    os.system("freebayes -f {} -F 0.15 -C 2 --pooled-continuous {}/{}_sorted.bam > {}/{}_freebayes.vcf".format(ref,out_file, sample,out_file, sample))
    #os.system("sed 's/ENA|MN908947|MN908947.3/NC_045512.2/g' {}/{}_freebayes.vcf > {}/{}_freebayes_sed.vcf".format(out_file, sample,out_file, sample))
    #os.system("java -jar snpEff/snpEff.jar ann NC_045512.2 {}/{}_freebayes_sed.vcf > {}/{}_spEff.vcf".format(out_file, sample,out_file, sample))
 
    
def fastq_pipeline(sample_dir, work_dir , ref, vcf = False):
    sample_number = os.listdir(sample_dir)[0].split("/")[-1].split("_")[0]    
    out_file = work_dir + "results/" + sample_number 
    os.system("mkdir {}".format(out_file))
    fs1out = out_file + "/" + sample_number + "R1_clean.fastq.gz"
    fs2out = out_file + "/" + sample_number + "R2_clean.fastq.gz"
    #os.system("fastp -i {} -I {} -o {} -O {} -h {}/{}_fastp.html -j {}/{}_fastp.json".format(sample_dir +"/"+ os.listdir(sample_dir)[0], sample_dir + "/" + os.listdir(sample_dir)[1], fs1out, fs2out,out_file, sample_number, out_file, sample_number))

    #os.system("bwa mem {} {} {} > {}/{}_aln_se.sam".format(ref, fs1out, fs2out, out_file, sample_number))
    ##### add a step where primers are removed
    #os.system("samtools sort {}/{}_aln_se.sam > {}/{}_sorted.bam".format(out_file, sample_number, out_file, sample_number))
    #os.system("samtools index {}/{}_sorted.bam".format(out_file, sample_number))
    #os.system("samtools consensus -f fasta {}/{}_sorted.bam -o {}/{}_consensus.fa ".format(out_file, sample_number, out_file, sample_number))
    #os.system("samtools depth -a {}/{}_sorted.bam -o {}/{}_coverage.txt".format(out_file, sample_number, out_file, sample_number)) 
    #os.system("samtools stats {}/{}_sorted.bam > {}/{}_stats.txt".format(out_file, sample_number, out_file, sample_number)) 
    
    if vcf :
        variant_calling(out_file,sample_number, ref)
    

    return "Sample {} processed!".format(sample_number) 
    
def coverage_plot(sample_dir):
    
    #res_dirs = os.listdir(sample_dir+"/results")
    #res_files = json_files = result_load(,"coverage.txt")
    res_files = result_load(sample_dir, "coverage")
    print(res_files)
    frames = []
    header = []
    empty = []
    for cv in res_files:
        if os.stat(cv).st_size != 0:  ##################################
            frames.append(pd.read_csv(cv, delimiter = "\t", usecols=[1,2], index_col = 0, names = ["position",cv.split("/")[-1].split("_")[0]]))
            header.append(cv.split("_")[0])
        else:
            empty.append(cv.split("_")[0])
            
    for e in empty:
        frames.append(pd.DataFrame(0, index = frames[0].index, columns = e))
        header.append(e)
    data = pd.concat(frames, axis = 1).sort_index(axis = 1)
    x = 0
    print(data)
    try:
        pdf = matplotlib.backends.backend_pdf.PdfPages("{}results/coverage_per_sample.pdf".format(sample_dir))
    except PermissionError:
        pdf = matplotlib.backends.backend_pdf.PdfPages("{}results/coverage output later {}.pdf".format(x))

    for i in range(0,data.shape[1]-4,8):
    #for i in range(0,8,8):    
        plt.rcParams.update({'font.size': 10})
    
        fig, ax = plt.subplots(2, 4, sharex='all', sharey='all' , figsize=(120, 60), constrained_layout=True)
        for q in range(0,4):
            ax[0,q].plot(data.index, np.log10(data.iloc[:,i+q]), color="#3399e6", lw=1 )
            ax[0,q].set_title("Coverage per base of sample No. {}".format(data.iloc[:,i+q].name), fontsize=60)
            ax[0,q].set_xlabel('Sequence Position', fontsize = 60)
            ax[0,q].set_ylabel('Coverage', fontsize = 60)
            ax[0,q].axhline(y=1, label='Quality Passed', c='black', lw = 2)
    
        for p in range(0,4):
            ax[1,p].plot(data.index, np.log10(data.iloc[:,i+4+p]), color="#3399e6", lw=1 )
            ax[1,p].set_title("Coverage per base of sample No. {}".format(data.iloc[:,i+4+p].name), fontsize=60)
            ax[1,p].set_xlabel('Sequence Position', fontsize = 60)
            ax[1,p].set_ylabel('Coverage', fontsize = 60)
            ax[1,p].axhline(y=1, label='Quality Passed', c='black', lw = 2)

    
        pdf.savefig( fig )

    pdf.close()
    
def result_load(input_file, file_type):
    res_dirs = os.listdir(input_file+"results")
    res_files = []
    for rd in res_dirs:
        if os.path.isdir(input_file+"results/" + rd):
            cd = os.listdir(input_file+"results/" + rd) 
            for cf in cd:
                if file_type in cf:   ######replace file type 
                    res_files.append(input_file+"results/" + rd +"/"+cf)
    return res_files
    
    
def quality_checks(file_folder):
    pd.set_option('display.float_format', lambda x: '%.4f' % x)
    data_tables = []
    json_files = result_load(file_folder,".json")
    for i in json_files:
        info = json.load(open(i))
        data_tables.append(pd.DataFrame(info["summary"]).add_prefix(i.split("/")[-1].split("_")[0] + "_"))
    return pd.concat(data_tables, axis = 1).sort_index()#.astype('int'))   

def read_distribution(file_folder):
    stat_files = result_load(file_folder,"stats")
    all_reads = {}
    for i in stat_files:
        sample = i.split("/")[-1].split("_")[0]
        f = ("/").join(i.split("/")[0:-1])
        os.system("grep ^RL {} > {}/{}_count_sizes.txt".format(i, f, sample))
        counts = pd.read_csv( "{}/{}_count_sizes.txt".format(f, sample), delimiter = "\t", names = ["Type","Size","Counts"] )
        all_reads[sample] = dict(zip(counts.Size, counts.Counts))
    return all_reads

def read_plot(file_folder):
    distributions = read_distribution(file_folder)
    try:
        pdf = matplotlib.backends.backend_pdf.PdfPages("{}results/read_length_per_sample.pdf".format(file_folder))
    except PermissionError:
        pdf = matplotlib.backends.backend_pdf.PdfPages("{}results/read_length_plains.pdf".format(file_folder))
    
    sample_list = sorted(list(distributions.keys()))
    for i in range(0,len(sample_list)-4,8):
    #for i in range(0,8,8):    
        plt.rcParams.update({'font.size': 10})
    
        fig, ax = plt.subplots(2, 4, sharex='all', sharey='all' , figsize=(120, 60),constrained_layout=True)
        for q in range(0,4):
            ax[0,q].bar(list(distributions[sample_list[i+q]].keys()), distributions[sample_list[i+q]].values(), color="#3399e6" , width = 0.6)
            ax[0,q].set_title("Read length distribution of sample No. {}".format(sample_list[i+q+4]), fontsize=60)
            ax[0,q].set_xlabel('Read length', fontsize = 60)
            ax[0,q].set_ylabel('Number of reads', fontsize = 60)
            #ax[0,q].axhline(y=1, label='Quality Passed', c='black', lw = 2)
    
        for p in range(0,4):
            ax[1,p].bar(list(distributions[sample_list[i+q+4]].keys()), distributions[sample_list[i+q+4]].values(), color="#3399e6" , width = 0.6)
            ax[1,p].set_title("Read length distribution of sample No. {}".format(sample_list[i+q+4]), fontsize=60)
            ax[1,p].set_xlabel('Read length', fontsize = 60)
            ax[1,p].set_ylabel('Number of reads', fontsize = 60)
            #ax[1,p].axhline(y=1, label='Quality Passed', c='black', lw = 2)

    
        pdf.savefig( fig )

    pdf.close()
        

def priming(h, bed):
    marker = 0
    for p in range(0,len(bed.index)):
        if h >= bed.iloc[p]["START"] and h <= bed.iloc[p]["END"]:
            return [0 ,"Hetero site found {}".format(h)]
    return [1,"Hetero site found {} in non priming region".format(h)]

def heterozygosity(sample_dir, primer_regions):
    het_data = {}
    os.system("mkdir {}/heterozygosity_stats".format(sample_dir))
    res_files = result_load(sample_dir, "freebayes.vcf")
    for i in res_files:
        f = open ( "{}/heterozygosity_stats/{}_heterozygosity.txt".format(sample_dir, i.split("/")[-1].split("_")[0]) , "w" )
        het_data[i.split("/")[-1].split("_")[0]] = 0
        callset = allel.read_vcf(i)
        gt = allel.GenotypeArray(callset['calldata/GT'])
        for g in range(0,len(gt)):
            #print(gt)
            if gt[g].is_het() and priming(callset['variants/POS'][g], primer_regions)[0] == 1 :
            #print(, g+1, callset['variants/POS'][g], callset['variants/REF'][g])
         
                f.write(i + "\t" + priming(callset['variants/POS'][g], primer_regions)[1] + "\t" + callset['variants/REF'][g] +"\n" )
                het_data[i.split("/")[-1].split("_")[0]]+=1
        f.close()  

    het_file = open(sample_dir +"results/heterozygosity_stats.csv", "w")
    for q in sorted(het_data.keys()):
        het_file.write(str(q) + "\t" + str(het_data[q]) + "\n")
    het_file.close()
    return het_data


################################################################

parser = argparse.ArgumentParser(description='A fast, efficient way to process paired-end data',  epilog = "author: Maria Malliarou <maria.malliarou.ger@gmail.com> v1.1" )

parser.add_argument('--input_path',  type = str, required = True, help = "Please provide your fastq_folder" )  
parser.add_argument('--primers',  type = str, help = "Please provide your bed file with the primer regions")  ###θα παίρνει ένα ή πολλα vcf αρχείο
parser.add_argument('--pools',  type = int, default = 2 ,help = "Please provide the number of multiprocess pools, if not the default is two") 
parser.add_argument('--reference',  type = str, default = "reference/GCA_009858895.3.fasta" ,help = "Please provide your reference fasta file, if not the program will use the default file") 
parser.add_argument('--vcf', default = False, type = bool, help = "If True the programm will also perform variant calling" )
parser.add_argument('--hetero_report', default = False, type = bool, help = "If True the programm will provide a report file with heterozygous sites in sites." )

args = parser.parse_args()
################################################################

if args.input_path.endswith("/") == False:
    args.input_path = args.input_path + "/"

files = [ fastq_files(args.input_path + f) for f in os.listdir(args.input_path) if f != "results" ]

items = [( f, args.input_path, args.reference , args.vcf) for f in files]

os.system("mkdir {}results".format(args.input_path))

#for i in items:
#    print(i)
#    fastq_pipeline(*i)


'''
with Pool(args.pools) as pool:
# call the same function with different data in parallel
    for result in pool.starmap(fastq_pipeline, items):
        print(result)



'''
if args.hetero_report:
    primer_bed = pd.read_csv("V4_primer_set.bed", delimiter = "\t", header = None, names = ["ID", "START", "END", "DIRECTION","DIR_ID","STRAND","SEQUENCE"])

    het_dict = heterozygosity(args.input_path, primer_bed)
    #het_file = open(args.input_path +"results/heterozygosity_stats.csv", "w")
    print(het_dict)
    #for sample in sorted(het_dict.keys()):
    #    het_file.write(sample + "\t" + het_dict[sample])
    #het_file.close()
    
    
#print(read_plot(args.input_path))
#coverage_plot(args.input_path)
#quality_checks(args.input_path).to_csv(args.input_path +"results/read_stats.csv")