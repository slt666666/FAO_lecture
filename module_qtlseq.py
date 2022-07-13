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

    mut_effect = [random.uniform(-2.5, 2.5) for i in range(len(mut_pos))]
    mut_effect[12] = 10

    pd.DataFrame({"Position": mut_pos, "Simulated_SNP_effect":mut_effect}).to_csv("../SNP_effect.csv")

    children = []
    phenotype = []
    for i in range(200):
        child = list(copy.copy(cultivar_A))
        phen = 40
        key = False
        for mut_i, j in enumerate(mut_pos):
            if random.random() > 0.8:
                key = not key
            if key:
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
    row = 1
    col = 2
    fig, ax=plt.subplots(row, col, figsize=(12,4))
    for i, SNP_index in enumerate([high_bulk_SNP_index, low_bulk_SNP_index]): 
        x = SNP_index.index.values
        y = SNP_index.SNP_index.values
        color = ["green", "orange"]; titles = ["High bulk SNP-index", "Low bulk SNP-index"]
        ax[i].scatter(x, y, color=color[i])
        ax[i].set_title(titles[i], fontsize=16)
        ax[i].set_ylim(-0.05, 1.05)
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

    tmp = delta_SNP_index[delta_SNP_index["delta_SNP_index"] >= 0.9]
    x = tmp.index.values
    y = tmp.delta_SNP_index.values
    plt.scatter(x, y, color="red")
    plt.ylim(-1.05, 1.05)
    plt.show()

def check_results(delta_SNP_index):
    SNP_effect = pd.read_csv("../SNP_effect.csv", index_col=0)
    SNP_effect = SNP_effect.set_index("Position")
    SNP_effect["delta_SNP_index"] = delta_SNP_index.delta_SNP_index.values
    return SNP_effect.reset_index()

def get_yesterday_SNP_index():
    !wget -q -O mutmap_dataset.txt https://raw.githubusercontent.com/CropEvol/lecture/master/data/mutmap_chr10.txt
    SNP_index = pd.read_csv("mutmap_dataset.txt", sep=',', header=0)
    SNP_index["SNP_index"] = SNP_index["N_ALT"] / (SNP_index["N_REF"] + SNP_index["N_ALT"])
    return SNP_index

def visualize_SNP_index2(SNP_index):
    x = SNP_index["POS"]
    y = SNP_index["SNP_index"]

    plt.figure(figsize=[12,4])
    plt.scatter(x, y)
    plt.title('SNP-index on chromosome 10', fontsize=18)
    plt.xlabel('Position (x 10 Mb)', fontsize=12)
    plt.ylabel('SNP-index', fontsize=12)

    df2 = SNP_index[ SNP_index["SNP_index"]>0.8 ]

    x2 = df2["POS"]
    y2 = df2["SNP_index"]
    plt.scatter(x2, y2, color="red")

    plt.show()   

def sliding_window(SNP_index, window_size=1 * 1000 * 1000, step_size = 0.2 * 1000 * 1000):
    chrom_size = 23207287

    x_win = []
    x = []
    y_win = []

    start = 1
    end  = start + window_size - 1

    while True:
        sub = SNP_index[ (SNP_index["POS"]>=start) & (SNP_index["POS"]<=end) ]
        p = sub["POS"].mean()
        s = sub["SNP_index"].mean()
        x_win.append(p)
        x.append([start, end])
        y_win.append(s)

        start = start + step_size
        end  = end  + step_size

        if end > chrom_size:
            break

    plt.figure(figsize=[12,4])
    plt.scatter(SNP_index["POS"], SNP_index["SNP_index"])
    plt.plot(x_win, y_win, color="red")
    plt.title('SNP-index on chromosome 10', fontsize=18)
    plt.xlabel('Position (x 10 Mb)', fontsize=12)
    plt.ylabel('SNP-index', fontsize=12)
    plt.show()

    result = pd.DataFrame({"Window":x, "average SNP-index":y_win})
    return result