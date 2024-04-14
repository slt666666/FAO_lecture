import os
import random
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
import matplotlib.patches as mpatches
from sklearn.linear_model import ElasticNet
sns.set()

def linear_model(X, y):
    model = ElasticNet(alpha=0.3, l1_ratio=0.7)
    model = model.fit(X, y["Grain_number"])
    equation = "Grain number = "
    k = 0
    for i, j in enumerate(model.coef_):
        if j == 0:
            pass
        else:
            k += 1
            equation += "SNP{}×{} + ".format(str(i+1), str(j))
    equation += str(model.intercept_)
    print(equation)
    return model

def check_accuracy(model, X, y):
    y_test_preds = model.predict(X)
    plt.figure(figsize=[4,4])
    plt.scatter(y.Grain_number, y_test_preds)
    plt.xlabel("Observed phenotype values")
    plt.ylabel("Predicted phenotype values")
    plt.title("correlation coefficient: {}".format(np.corrcoef(y.Grain_number, y_test_preds)[0, 1]), fontsize="large")
    plt.show()

def predict_phenotype(test_genotype, prediction_model):
    y_test_pred = prediction_model.predict(test_genotype.T.fillna(1))
    return y_test_pred

def predict_progeny_phenotype(Line1, Line2, progeny, phenotype, genotype, prediction_model):

    Line1_pheno = phenotype.loc[phenotype.Line == Line1, "Grain_number"].values[0]
    Line2_pheno = phenotype.loc[phenotype.Line == Line2, "Grain_number"].values[0]
    print(f"Grain number of {Line1} is {Line1_pheno}")
    print(f"Grain number of {Line2} is {Line2_pheno}")

    print("If we cross {} & {}, the phenotype of F2 population may be...".format(Line1, Line2))

    progenies = []

    for k in range(progeny):
        a = np.arange(len(genotype[Line1].values))
        a = np.random.permutation(a)[:len(genotype[Line1].values)//100]
        a.sort()

        if np.random.rand() > 0.5:
            new = copy.copy(genotype[Line1].values)
            other = copy.copy(genotype[Line2].values)
        else:
            new = copy.copy(genotype[Line2].values)
            other = copy.copy(genotype[Line1].values)

        for i in range(len(a)//2):
            new[a[i*2]:a[i*2+1]] = other[a[i*2]:a[i*2+1]]
        
        new[new < 0] = 0
        progenies.append(new)
    progenies_genotype = pd.DataFrame(progenies).T

    pred = predict_phenotype(progenies_genotype, prediction_model)
    sns.set()
    plt.figure(figsize=[4,4])
    plt.hist(pred)
    plt.show()

def predict_customized_genotype(genotype, regions, prediction_model):
    trait = "GN"
    fig = plt.figure(figsize=(5,12))
    ax = plt.axes()

    for i, each_chr in enumerate(genotype.chr.unique()):
        chr_genotype = genotype[genotype["chr"] == each_chr]
        end = chr_genotype.iloc[-1, :].pos
        r = patches.Rectangle(xy=(i*16000000, 0), width=6000000, height=-end*3, ec='gray', fc="orange", linewidth=3)
        ax.add_patch(r)
        if i == 0:
            plt.text(-10000000, 2000000, "chr {}".format(i+1))
        else:
            plt.text(i*16000000, 2000000, str(i+1))

    r = patches.Rectangle(xy=(65000000, -1*130000000), width=3000000*3, height=5000000, ec='gray', fc="orange", linewidth=3)
    ax.add_patch(r)
    plt.text(29000000, -1*130000000, "Cutivar A")
    r = patches.Rectangle(xy=(142000000, -1*130000000), width=3000000*3, height=5000000, ec='gray', fc="blue", linewidth=3)
    ax.add_patch(r)
    plt.text(80000000, -1*130000000, "Another Cultivar")

    for region in regions:
        i = int(region[0][3:]) - 1
        chr_genotype = genotype[genotype["chr"] == region[0]]
        end = chr_genotype.iloc[-1, :].pos
        if region[1] >= region[2]:
            print("region is uncorrect.")
            break
        elif region[1] < 0:
            print("start sould be larger than 0.")
            break
        elif region[2] > end:
            regions[2] = end
        else:
            r = patches.Rectangle(xy=(i*16000000, -3*region[1]), width=6000000, height=-(region[2] - region[1])*3, ec='gray', fc="blue", linewidth=3)
            ax.add_patch(r)

    plt.axis('scaled')
    plt.axis('off')
    plt.tight_layout()
    ax.set_aspect('equal')

    plt.show()

    customized_genotype = np.repeat(0, genotype.shape[0])
    for region in regions:
        selected_chr = region[0]
        selected_start = region[1]
        selected_end = region[2]
        customized_genotype[genotype.chr.isin([selected_chr]) & (genotype.pos >= selected_start) & (genotype.pos <= selected_end)] = 2
    
    print("\n")
    print("The cultivar of this new genotype showed...")
    print(trait, "=", '\033[1m'+"{}".format(predict_phenotype(pd.DataFrame(customized_genotype), prediction_model)[0])+'\033[0m')
    print("※", trait, "of Cultivar A is", '\033[1m'+"{}".format(predict_phenotype(pd.DataFrame(np.repeat(0, genotype.shape[0])), prediction_model)[0])+'\033[0m')