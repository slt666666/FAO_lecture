import os
import random
import pandas as pd
import numpy as np
import igv_notebook
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def RefTrack(track_info, locus=None):
    ref ={"reference": {**{"id":"Reference"}, **track_info}, "locus": locus}
    return ref

def BamTrack(track_info):
    base = {"name":None, "type":"alignment", "format":"bam"}
    bam = {**base, **track_info}
    return bam

def AnnotationTrack(track_info):
    base = {"name":None, "type":"annotation", "format":"gff"}
    bam = {**base, **track_info}
    return bam

def calc_SNP_index(vcf):
    vcf = pd.read_csv(vcf, comment="#", sep="\t", header=None)
    vcf_data = []
    for i in range(vcf.shape[0]):
        position = vcf.iloc[i, 1]
        reference =  vcf.iloc[i, 3]
        mutation =  vcf.iloc[i, 4]
        read_num = vcf.iloc[i, 9][vcf.iloc[i, 9].rfind(":")+1:].split(",")
        ref_reads = int(read_num[0])
        mut_reads = int(read_num[1])
        total_reads = ref_reads+mut_reads
        alt_per = mut_reads / total_reads
        vcf_data.append([position, reference, mutation, mut_reads, total_reads, alt_per])
    vcf_data = pd.DataFrame(vcf_data)
    vcf_data.columns = ["position", "ref", "mutant", "mutant base", "total reads", "SNP index"]
    return vcf_data

def calc_delta_SNP_index(vcf):
    vcf = pd.read_csv(vcf, comment="#", sep="\t", header=None)
    vcf_data = []
    for i in range(vcf.shape[0]):
        position = vcf.iloc[i, 1]
        reference =  vcf.iloc[i, 3]
        mutation =  vcf.iloc[i, 4]
        low_read_num = vcf.iloc[i, 9][vcf.iloc[i, 9].rfind(":")+1:].split(",")
        low_ref_reads = int(low_read_num[0])
        low_mut_reads = int(low_read_num[1])
        low_total_reads = low_ref_reads+low_mut_reads
        low_SNP_index = low_mut_reads / low_total_reads
        high_read_num = vcf.iloc[i, 10][vcf.iloc[i, 10].rfind(":")+1:].split(",")
        high_ref_reads = int(high_read_num[0])
        high_mut_reads = int(high_read_num[1])
        high_total_reads = high_ref_reads+high_mut_reads
        high_SNP_index = high_mut_reads / high_total_reads
        delta = high_SNP_index - low_SNP_index
        vcf_data.append([position, reference, mutation, low_SNP_index, high_SNP_index, delta])
    vcf_data = pd.DataFrame(vcf_data)
    vcf_data.columns = ["position", "ref", "mutant", "low SNP index", "high SNP index", "delta SNP index"]
    return vcf_data

def visualize_SNP_index(vcf_data):
    plt.figure(figsize=(6, 2))
    sns.scatterplot(x=vcf_data.position, y=vcf_data["SNP index"], color="black")
    plt.axhline(y=0.5, color="red")
    plt.ylim(-0.1, 1.1)
    plt.show()

def visualize_delta_SNP_index(vcf_data):
    plt.figure(figsize=(6, 2))
    sns.scatterplot(x=vcf_data.position, y=vcf_data["low SNP index"], color="orange")
    plt.axhline(y=0.5, color="red")
    plt.ylim(-0.1, 1.1)
    plt.xlabel("")
    plt.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
    plt.show()
    plt.figure(figsize=(6, 2))
    sns.scatterplot(x=vcf_data.position, y=vcf_data["high SNP index"], color="green")
    plt.axhline(y=0.5, color="red")
    plt.ylim(-0.1, 1.1)
    plt.xlabel("")
    plt.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
    plt.show()
    plt.figure(figsize=(6, 2))
    sns.scatterplot(x=vcf_data.position, y=vcf_data["delta SNP index"], color="black")
    plt.axhline(y=0, color="red")
    plt.ylim(-1.1, 1.1)
    plt.show()

def simulate_fastq(F2_num, read_num):
    with open("genome/CultivarB.fa", mode="r") as f:
        reference = f.readlines()[1]
    with open("simulation/mutations.fa", mode="r") as f:
        mutation = f.readlines()[1]

    F2_genotype = [list(reference)]
    F2_seq = []
    for i in range(F2_num):
        recom_points = np.sort([random.randint(0, 9999) for i in range(10)])
        if random.randint(0, 9999) > 5000:
            new_geno = reference[:recom_points[0]]+mutation[recom_points[0]:recom_points[1]]+reference[recom_points[1]:recom_points[2]]+mutation[recom_points[2]:recom_points[3]]+\
            reference[recom_points[3]:recom_points[4]]+mutation[recom_points[4]:recom_points[5]]+reference[recom_points[5]:recom_points[6]]+mutation[recom_points[6]:recom_points[7]]+\
            reference[recom_points[7]:recom_points[8]]+mutation[recom_points[8]:recom_points[9]]+reference[recom_points[9]:]
        else:
            new_geno = mutation[:recom_points[0]]+reference[recom_points[0]:recom_points[1]]+mutation[recom_points[1]:recom_points[2]]+reference[recom_points[2]:recom_points[3]]+\
            mutation[recom_points[3]:recom_points[4]]+reference[recom_points[4]:recom_points[5]]+mutation[recom_points[5]:recom_points[6]]+reference[recom_points[6]:recom_points[7]]+\
            mutation[recom_points[7]:recom_points[8]]+reference[recom_points[8]:recom_points[9]]+mutation[recom_points[9]:]
        F2_genotype.append(list(new_geno))
        F2_seq.append(new_geno)

    F2_genotype = pd.DataFrame(F2_genotype)
    F2_genotype = F2_genotype.loc[:, (F2_genotype == F2_genotype.iloc[0, :]).sum() != F2_num+1]

    with open("simulation/F2_genome.fa", "w") as f:
        for i, key in enumerate((F2_genotype.iloc[1:, 19] != F2_genotype.iloc[0, 19]).values):
            if key:
                f.write(f">sample{i}\n")
                f.write(F2_seq[i])
                f.write("\n")

    os.system(f'wgsim -e 0 -r 0 -R 0 -X 0 -d 300 -1 150 -2 150 -N {read_num} simulation/F2_genome.fa simulation/bulked_1.fastq simulation/bulked_2.fastq')
