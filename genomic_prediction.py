import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNet

def load_dataset():
    phenotype = pd.read_csv("")
    genotype = pd.read_csv("")
    return genotype, phenotype

def split_dataset(genotype, phenotype, trait, test=0.2):
    train, test = train_test_split(phenotype.Line.values, test_size=0.2)
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

def predict_phenotype(test_genotype, prediction_model):
    y_test_pred = clf.predict(test_genotype.T.fillna(1))
    return y_test_pred

def check_accuracy(predicted_test_phenotype, test_data, trait):
    
    print("Correlation coefficient:", np.corrcoef([y_test_pred, test_phenotype[trait].values])[0,1])
    sns.set()
    plt.scatter(y_test_pred, test_phenotype[trait].values)
    plt.xlabel("predicted {}".format(trait))
    plt.ylabel("observed {}".format(trait))
    plt.show()