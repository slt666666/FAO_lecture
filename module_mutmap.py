import random
import copy
import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot  as plt
import seaborn as sns
import requests
import io

def make_mutant(reference):
    mutant = list(copy.copy(reference))
    mut_pos = [5]
    l = list(range(len(reference)))
    l.remove(5)
    mut_pos.extend(random.sample(l, 4))
    for i in range(5):
        mutants = ["A", "T", "G", "C"]
        mutants.remove(reference[mut_pos[i]])
        mutant[mut_pos[i]] = random.choice(mutants)
    mutant = "".join(mutant)

    print("reference genotype: ", reference)
    print("   mutant genotype: ", mutant)

    with open("reference_sequences.fasta", mode='w') as f:
        f.write(">reference\n")
        f.write("{}\n".format(reference))
    
    return mutant

def cross_reference_and_mutant(reference, mutant, progeny=200):
    children = []
    phenotype = []
    for i in range(progeny):
        child = ""
        for j in range(len(reference)):
            if random.random() > 0.5:
                child = child + reference[j]
            else:
                child = child + mutant[j]
        children.append(child)
        if child[5] == reference[5]:
            phenotype.append("Green")
        else:
            phenotype.append("LightGreen")
    children = pd.DataFrame(children)
    children.columns = ["Genotype"]
    children["Phenotype"] = phenotype
    children.to_csv("../simulated_progenies.csv")

    image_ref = Image.open(io.BytesIO(requests.get('https://github.com/slt666666/FAO_lecture/blob/main/reference_plant.png?raw=true').content))
    image_mut = Image.open(io.BytesIO(requests.get('https://github.com/slt666666/FAO_lecture/blob/main/mutant_plant.png?raw=true').content))
    image_dash = Image.open(io.BytesIO(requests.get('https://github.com/slt666666/FAO_lecture/blob/main/dash.png?raw=true').content))

    images = []
    for i in range(children.shape[0]):
        if children.iloc[i, 1] == "LightGreen":
            images.append(image_mut)
        else:
            images.append(image_ref)
    print("Progenies are ...")
    fig = plt.figure()
    for i, im in enumerate(images[:20]):
        if i != 19:
            fig.add_subplot(4,5,i+1).set_title(str(i))
            plt.axis("off")
            plt.imshow(im)
        else:
            fig.add_subplot(4,5,i+1).set_title(str(i)+"...")
            plt.axis("off")
            plt.imshow(image_dash)
    plt.show()

    return children

def bulk_sequencing(progeny):
    mutants = progeny[progeny["Phenotype"] == "LightGreen"]
    mutants.iloc[:, 0].values

    reads = []
    for i in range(50):
        tmp = np.random.choice(mutants.iloc[:, 0].values)
        start = np.random.choice(range(int(len(tmp)//2)+1))
        read = tmp[start:start+10]
        reads.append([read, start])

    reads = pd.DataFrame(reads)
    reads.to_csv("../reads.csv")

    with open("bulked_sequences.fasta", mode='w') as f:
        for i in range(reads.shape[0]):
            f.write(">read{}\n".format(i))
            f.write("{}\n".format(reads.iloc[i, 0]))
    
    return reads

def alignment(reads, reference):    
    tmp_df = pd.DataFrame(list(reference))
    tmp_df.columns = ["Reference"]
    reads = reads.sort_values(by=1)
    for i in range(reads.shape[0]):
        tmp_df["read{}".format(i)] = list("-"*reads.iloc[i, 1] + reads.iloc[i, 0] + "-"*(10-reads.iloc[i, 1]))
    tmp_df = tmp_df.T
    display(tmp_df)
    return tmp_df

def calculate_SNP_index(alignment_result, reference, mutant):
    table_df = pd.DataFrame(list(reference))
    table_df.columns = ["Reference"]
    table_df["Mutant"] = pd.DataFrame(list(mutant))
    ref_num = []
    mut_num = []
    for i in range(alignment_result.shape[1]):
        tmp = alignment_result.iloc[1:, i]
        tmp = list(tmp[tmp != "-"].values)
        if reference[i] != mutant[i]:
            ref_num.append(tmp.count(reference[i]))
            mut_num.append(tmp.count(mutant[i]))
        else:
            ref_num.append(tmp.count(reference[i]))
            mut_num.append(0)
    table_df["ref_num"] = ref_num
    table_df["mut_num"] = mut_num
    table_df["SNP_index"] = table_df["mut_num"] / (table_df["mut_num"] + table_df["ref_num"])
    display(table_df)
    return table_df

def visualize_SNP_index(SNP_index):
    sns.set()
    x = SNP_index.index.values
    y = SNP_index.SNP_index.values
    plt.xticks(range(SNP_index.shape[0]), range(1, SNP_index.shape[0]+1))
    plt.scatter(x, y)
    plt.show()

def check_results(reference):
    reference_df = pd.DataFrame(list(reference)).T
    reference_df.columns = range(1, 21)
    reference_df["Phenotype"] = "Green"
    reference_df.index = ["Reference"]
    display(reference_df)

    children = pd.read_csv("../simulated_progenies.csv", index_col=0)
    children_sample1 = pd.DataFrame(children[children["Phenotype"] == "LightGreen"].iloc[:5, 0].str.split("").tolist())
    children_sample1 = children_sample1.iloc[:, 1:-1]
    children_sample1["Phenotype"] = "LightGreen"
    children_sample2 = pd.DataFrame(children[children["Phenotype"] == "Green"].iloc[:5, 0].str.split("").tolist())
    children_sample2 = children_sample2.iloc[:, 1:-1]
    children_sample2["Phenotype"] = "Green"
    children_sample = pd.concat([children_sample1, children_sample2])
    children_sample.index = ["children"+str(i) for i in range(children_sample.shape[0])]
    display(children_sample)