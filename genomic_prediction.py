import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNet

def load_dataset():
    phenotype = pd.read_csv("https://github.com/slt666666/FAO_lecture/blob/main/phenotype.csv?raw=true")
    genotype = pd.read_csv("https://github.com/slt666666/FAO_lecture/blob/main/genotype.csv?raw=true")
    return genotype, phenotype

def split_dataset(genotype, phenotype, trait, test=0.2):
    train, test = train_test_split(phenotype.Line.values, test_size=0.2, random_state=1)
    train = np.sort(train); test = np.sort(test)
    test_genotype = genotype.loc[:, test]
    test_phenotype = phenotype[phenotype.Line.isin(test)]
    train_genotype = genotype.loc[:, train]
    train_phenotype = phenotype[phenotype.Line.isin(train)]
    return test_genotype, test_phenotype, train_genotype, train_phenotype

def make_genomic_prediction_model(train_genotype, train_phenotype, phenotype):
    
    X_train = train_genotype.T.fillna(1)
    y_train = train_phenotype[phenotype].values
    if phenotype == "LW_mean":
        alpha=0.02
        l1_ratio=0.1
    else:
        alpha=0.5
        l1_ratio=0.3
    clf = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
    clf.fit(X_train, y_train)
    # coefs = clf.coef_
    return clf

def check_equation(trait, prediction_model):
    equation = "{} = {} + ".format(trait, prediction_model.intercept_)
    for i, j in enumerate(prediction_model.coef_):
        if j != 0:
            equation += "SNP{}×{} + ".format(str(i), str(j))
    return equation

def predict_phenotype(test_genotype, prediction_model):
    y_test_pred = prediction_model.predict(test_genotype.T.fillna(1))
    return y_test_pred

def check_accuracy(predicted_test_phenotype, test_phenotype, trait):
    
    print("Correlation coefficient:", np.corrcoef([predicted_test_phenotype, test_phenotype[trait].values])[0,1])
    sns.set()
    plt.scatter(predicted_test_phenotype, test_phenotype[trait].values)
    plt.xlabel("predicted {}".format(trait))
    plt.ylabel("observed {}".format(trait))
    plt.show()

def show_estimated_SNP_effect(prediction_model):
    sns.set()
    plt.plot(prediction_model.coef_)
    plt.show()

def predict_progeny_phenotype(Line1, Line2, progeny, genotype, prediction_model):

    print("If we cross {} & {}, the phenotype of F2 population may be...".format(Line1, Line2))

    genotype[Line1].values
    naguruu = genotype[Line2].values

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

        progenies.append(new)
    progenies_genotype = pd.DataFrame(progenies).T

    pred = predict_phenotype(progenies_genotype, prediction_model)
    sns.set()
    plt.hist(pred)
    plt.show()