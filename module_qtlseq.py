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