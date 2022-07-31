import random
import copy
import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot  as plt
import seaborn as sns
import requests
import io

def make_reference_and_mutant(length=20, mutation=5):
    reference = "".join([random.choice(["A", "T", "G", "C"]) for i in range(length)])
    mutant = list(copy.copy(reference))

    l = list(range(4, len(reference)-3))
    mut_pos = random.sample(l, mutation)
    for i in range(mutation):
        mutants = ["A", "T", "G", "C"]
        mutants.remove(reference[mut_pos[i]])
        mutant[mut_pos[i]] = random.choice(mutants)
    mutant = "".join(mutant)

    with open("reference_sequences.fasta", mode='w') as f:
        f.write(">reference\n")
        f.write("{}\n".format(reference))
    with open("../mutant_sequences.fasta", mode='w') as f:
        f.write(">mutant\n")
        f.write("{}\n".format(mutant))
    
    return reference, mutant

def cross_reference_and_mutant(reference, mutant, progeny=200):
    mut_pos = []
    for i, j in enumerate(list(reference)):
        if j != mutant[i]:
            mut_pos.append(i)
    causative_pos = random.choice(mut_pos)

    children = []
    phenotype = []
    for i in range(progeny):
        child = ""
        key = False
        for j in range(len(reference)):
            if random.random() > 0.8:
                key = not key
            if key:
                child = child + mutant[j]
            else:
                child = child + reference[j]
        children.append(child)
        if child[causative_pos] == reference[causative_pos]:
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
    print("F2 population is ...")
    fig = plt.figure(figsize=(24, 20))
    if progeny > 200:
        for i, im in enumerate(images[:200]):
            if i != 199:
                fig.add_subplot(10,20,i+1).set_title(str(i))
                plt.axis("off")
                plt.imshow(im)
            else:
                fig.add_subplot(10,20,i+1).set_title(str(i)+"...")
                plt.axis("off")
                plt.imshow(image_dash)
    elif progeny == 200:
        for i, im in enumerate(images[:200]):
            fig.add_subplot(10,20,i+1).set_title(str(i))
            plt.axis("off")
            plt.imshow(im)
    else:
        for i, im in enumerate(images[:20]):
            if i != 19:
                fig.add_subplot(4,5,i+1).set_title(str(i+1))
                plt.axis("off")
                plt.imshow(im)
            else:
                fig.add_subplot(4,5,i+1).set_title(str(i+1)+"...")
                plt.axis("off")
                plt.imshow(image_dash) 
    plt.show()

    return children

def bulk_sequencing(progeny, read=50):
    mutants = progeny[progeny["Phenotype"] == "LightGreen"]
    mutants.iloc[:, 0].values

    reads = []
    for i in range(read):
        tmp = np.random.choice(mutants.iloc[:, 0].values)
        start = np.random.choice(range(len(tmp)-19))
        read = tmp[start:start+20]
        reads.append([read, start])

    reads = pd.DataFrame(reads)
    reads.to_csv("../reads.csv")

    with open("bulked_sequences.fastq", mode='w') as f:
        for i in range(reads.shape[0]):
            f.write(">read{}\n".format(i))
            f.write("{}\n".format(reads.iloc[i, 0]))
            f.write("+\n")
            quality = random.choices(['A', 'B', 'J', 'I', '?', '&', '9', '7'], k=8)
            f.write("?J"+"".join(quality)+"\n")
    
    return reads

def alignment(reads, reference):    
    tmp_df = pd.DataFrame(list(reference))
    tmp_df.columns = ["Reference"]
    reads = reads.sort_values(by=1)
    for i in range(reads.shape[0]):
        tmp_df["read{}".format(reads.index.values[i])] = list("-"*reads.iloc[i, 1] + reads.iloc[i, 0] + "-"*(len(reference)-20-reads.iloc[i, 1]))
    tmp_df = tmp_df.T
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
    table_df.index = table_df.index.values + 1
    table_df = table_df[table_df["mut_num"] > 0]
    return table_df

def visualize_SNP_index(SNP_index):
    sns.set()
    x = SNP_index.index.values
    y = SNP_index.SNP_index.values
    plt.scatter(x, y)
    plt.xticks(np.arange(5, SNP_index.index.values[-1], 5))
    plt.xlabel("SNP position")
    plt.ylabel("SNP-index")
    plt.show()

def check_results(reference, SNP_index):
    reference_df = pd.DataFrame(list(reference)).T
    reference_df.columns = range(1, len(reference)+1)
    reference_df["Phenotype"] = "Green"
    reference_df.index = ["Reference"]

    children = pd.read_csv("../simulated_progenies.csv", index_col=0)
    children_sample1 = pd.DataFrame(children[children["Phenotype"] == "LightGreen"].iloc[:5, 0].str.split("").tolist())
    children_sample1 = children_sample1.iloc[:, 1:-1]
    children_sample1["Phenotype"] = "LightGreen"
    children_sample2 = pd.DataFrame(children[children["Phenotype"] == "Green"].iloc[:5, 0].str.split("").tolist())
    children_sample2 = children_sample2.iloc[:, 1:-1]
    children_sample2["Phenotype"] = "Green"
    children_sample = pd.concat([children_sample1, children_sample2])
    children_sample.index = ["children"+str(i) for i in range(children_sample.shape[0])]
    
    final_df = pd.concat([reference_df, children_sample])
    return final_df#.loc[:, [SNP_index.index[SNP_index["SNP_index"] == 1], "Phenotype"]]

def load_data():
    dataset = "mutmap_dataset.txt"
    df = pd.read_csv(dataset)
    return df

def calculate_SNP_index2(alignment_results):
    alignment_results["SNP_index"] = alignment_results["N_ALT"] / (alignment_results["N_REF"] + alignment_results["N_ALT"])
    return alignment_results

def visualize_SNP_index2(SNP_index):
    sns.set()
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

def MutMap_simulation(length=100, mutation=20, progeny=200, read=1000):
    reference, mutant = make_reference_and_mutant(length=length, mutation=mutation)
    progeny = cross_reference_and_mutant(reference, mutant, progeny=progeny)
    reads = bulk_sequencing(progeny, read=read)
    alignment_result = alignment(reads, reference)
    SNP_index = calculate_SNP_index(alignment_result, reference, mutant)
    display(SNP_index)
    visualize_SNP_index(SNP_index)
    check_results(reference)