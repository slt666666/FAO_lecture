import pandas as pd
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