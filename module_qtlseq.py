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