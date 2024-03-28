import pandas as pd
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

def visualize_SNP_index(vcf_data):
    plt.figure(figsize=(6, 2))
    sns.scatterplot(x=vcf_data.position, y=vcf_data["SNP index"], color="black")
    plt.axhline(y=0.5, color="red")
    plt.ylim(-0.1, 1.1)
    plt.show()