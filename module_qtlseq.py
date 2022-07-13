import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import copy

def make_2_cultivars(length=100, snp=20):
    cultivar_A = [random.choice(["A", "T", "G", "C"]) for i in range(length)]
    cultivar_A = "".join(cultivar_A)

    cultivar_B = list(copy.copy(cultivar_A))
    mut_pos = []
    l = list(range(len(cultivar_A)))
    mut_pos.extend(random.sample(l, snp))
    for i in range(snp):
        mutants = ["A", "T", "G", "C"]
        mutants.remove(cultivar_A[mut_pos[i]])
        cultivar_B[mut_pos[i]] = random.choice(mutants)
    cultivar_B = "".join(cultivar_B)
    return cultivar_A, cultivar_B

def make_F2_progeny(cultivar_A, cultivar_B, progeny=200):
    mut_pos = []
    for i, j in enumerate(list(cultivar_A)):
        if j != cultivar_B[i]:
            mut_pos.append(i)

    mut_effect = [random.uniform(-2.5, 2.5) for i in range(20)]
    mut_effect[5] = 10

    children = []
    phenotype = []
    for i in range(200):
        child = list(copy.copy(cultivar_A))
        phen = 40
        for mut_i, j in enumerate(mut_pos):
            if random.random() > 0.5:
                child[j] = cultivar_B[j]
                phen = phen + mut_effect[mut_i]
        phen += random.uniform(-1, 1)
        children.append("".join(child))
        phenotype.append(phen)

    children = pd.DataFrame(children)
    children.columns = ["Genotype"]
    children["Phenotype"] = phenotype
    return children

def check_distribution(f2_progeny):
    sns.set()
    plt.hist(f2_progeny["Phenotype"])

def high_and_low_bulk_sequencing(progeny, top=20, bottom=20, reads=500):
    children_sort = progeny.sort_values(by="Phenotype")
    high_bulk = children_sort.iloc[-1*top:, :]
    low_bulk = children_sort.iloc[:bottom, :]

    low_reads = []
    for i in range(reads):
        tmp = np.random.choice(low_bulk.iloc[:, 0].values)
        start = np.random.choice(range(len(tmp)-9))
        read = tmp[start:start+10]
        low_reads.append([read, start])

    low_reads = pd.DataFrame(low_reads)
    low_reads.to_csv("../low_reads.csv")

    with open("low_bulked_sequences.fasta", mode='w') as f:
        for i in range(low_reads.shape[0]):
            f.write(">read{}\n".format(i))
            f.write("{}\n".format(low_reads.iloc[i, 0]))

    high_reads = []
    for i in range(reads):
        tmp = np.random.choice(high_bulk.iloc[:, 0].values)
        start = np.random.choice(range(len(tmp)-9))
        read = tmp[start:start+10]
        high_reads.append([read, start])

    high_reads = pd.DataFrame(high_reads)
    high_reads.to_csv("../high_reads.csv")

    with open("high_bulked_sequences.fasta", mode='w') as f:
        for i in range(high_reads.shape[0]):
            f.write(">read{}\n".format(i))
            f.write("{}\n".format(high_reads.iloc[i, 0]))

    return high_reads, low_reads

def alignment(high_reads, low_reads, cultivar_A):
    low_df = pd.DataFrame(list(cultivar_A))
    low_df.columns = ["Reference"]
    low_reads = low_reads.sort_values(by=1)
    for i in range(low_reads.shape[0]):
        low_df["read{}".format(i)] = list("-"*low_reads.iloc[i, 1] + low_reads.iloc[i, 0] + "-"*(90-low_reads.iloc[i, 1]))
    low_df = low_df.T

    high_df = pd.DataFrame(list(cultivar_A))
    high_df.columns = ["Reference"]
    high_reads = high_reads.sort_values(by=1)
    for i in range(high_reads.shape[0]):
        high_df["read{}".format(i)] = list("-"*high_reads.iloc[i, 1] + high_reads.iloc[i, 0] + "-"*(90-high_reads.iloc[i, 1]))
    high_df = high_df.T

    return high_df, low_df

def calculate_SNP_index(high_bulk_alignment_result, low_bulk_alignment_result, cultivar_A, cultivar_B):
    low_table_df = pd.DataFrame(list(cultivar_A))
    low_table_df.columns = ["cultivarA"]
    low_table_df["cultivarB"] = pd.DataFrame(list(cultivar_B))
    ref_num = []
    mut_num = []
    for i in range(low_bulk_alignment_result.shape[1]):
        tmp = low_bulk_alignment_result.iloc[1:, i]
        tmp = list(tmp[tmp != "-"].values)
        if cultivar_A[i] != cultivar_B[i]:
            ref_num.append(tmp.count(cultivar_A[i]))
            mut_num.append(tmp.count(cultivar_B[i]))
        else:
            ref_num.append(tmp.count(cultivar_A[i]))
            mut_num.append(0)
    low_table_df["ref_num"] = ref_num
    low_table_df["mut_num"] = mut_num
    low_table_df["SNP_index"] = low_table_df["mut_num"] / (low_table_df["mut_num"] + low_table_df["ref_num"])

    high_table_df = pd.DataFrame(list(cultivar_A))
    high_table_df.columns = ["cultivarA"]
    high_table_df["cultivarB"] = pd.DataFrame(list(cultivar_B))
    ref_num = []
    mut_num = []
    for i in range(high_bulk_alignment_result.shape[1]):
        tmp = high_bulk_alignment_result.iloc[1:, i]
        tmp = list(tmp[tmp != "-"].values)
        if cultivar_A[i] != cultivar_B[i]:
            ref_num.append(tmp.count(cultivar_A[i]))
            mut_num.append(tmp.count(cultivar_B[i]))
        else:
            ref_num.append(tmp.count(cultivar_A[i]))
            mut_num.append(0)
    high_table_df["ref_num"] = ref_num
    high_table_df["mut_num"] = mut_num
    high_table_df["SNP_index"] = high_table_df["mut_num"] / (high_table_df["mut_num"] + high_table_df["ref_num"])
    
    high_table_df = high_table_df[high_table_df["cultivarA"] != high_table_df["cultivarB"]]
    low_table_df = low_table_df[low_table_df["cultivarA"] != low_table_df["cultivarB"]]
    return high_table_df, low_table_df

def visualize_SNP_index(high_bulk_SNP_index, low_bulk_SNP_index):
    sns.set()
    x = high_bulk_SNP_index.index.values
    y = high_bulk_SNP_index.SNP_index.values
    y2 = low_bulk_SNP_index.SNP_index.values
    plt.scatter(x, y)
    plt.title("High bulk SNP-index")
    plt.show()
    plt.clf()

    plt.scatter(x, y2)
    plt.title("Low bulk SNP-index")
    plt.show()


def calculate_delta_SNP_index(high_bulk_SNP_index, low_bulk_SNP_index):
    high_bulk_SNP_index["delta_SNP_index"] = high_bulk_SNP_index.SNP_index.values - low_bulk_SNP_index.SNP_index.values
    delta_SNP_index = high_bulk_SNP_index.loc[:, ["cultivarA", "cultivarB", "delta_SNP_index"]]
    return delta_SNP_index

def visualize_delta_SNP_index(delta_SNP_index):
    sns.set()
    x = delta_SNP_index.index.values
    y = delta_SNP_index.delta_SNP_index.values
    plt.scatter(x, y)
    plt.show()