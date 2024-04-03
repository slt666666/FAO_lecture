import os
import random
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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