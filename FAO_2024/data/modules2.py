import os
import random
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
    plt.scatter(y.Grain_number, y_test_preds)
    plt.xlabel("Observed phenotype values")
    plt.ylabel("Predicted phenotype values")
    plt.text(1, 0.5, "correlation coefficient: {}".format(np.corrcoef(y.Grain_number, y_test_preds)[0, 1]), fontsize="xx-large")