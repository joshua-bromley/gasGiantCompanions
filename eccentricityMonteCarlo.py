import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats

def eccenRedraw(eccens, errors):
    newEccens = np.zeros_like(eccens)
    for i in range(len(eccens)):
        newEccens[i] = stats.norm.rvs(loc = eccens[i], scale = errors[i])
    return newEccens

planets = pd.read_csv("./data/gasGiantData.csv")

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
seCompanions = planets.loc[planets["companion_type"] % 2 == 0]

gasGiantEccen = np.array(gasGiants["pl_orbeccen"])
gasGiantEccenErr = (np.array(gasGiants["pl_orbeccenerr1"]) - np.array(gasGiants["pl_orbeccenerr2"]))/2

seCompanionEccen = np.array(seCompanions["pl_orbeccen"])
seCompanionEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

ksResults = np.ones(50)

for i in range(len(ksResults)):
    gasGiantEccenRedraw = eccenRedraw(gasGiantEccen, gasGiantEccenErr)
    seCompanionEccenRedraw = eccenRedraw(seCompanionEccen, seCompanionEccenErr)
    ks = stats.kstest(seCompanionEccenRedraw, gasGiantEccenRedraw)
    ksResults[i] = ks[1]

print(np.mean(ksResults))
