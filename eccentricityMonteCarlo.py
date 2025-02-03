import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats
import corner

def eccenRedraw(eccens, errors):
    newEccens = np.zeros_like(eccens)
    for i in range(len(eccens)):
        if not (np.isnan(eccens[i]) or np.isnan(errors[i])):
            newEccens[i] = stats.norm.rvs(loc = eccens[i], scale = errors[i])
        elif np.isnan(errors[i]):
            while newEccens[i] <= 0:
                newEccens[i] = stats.norm.rvs(loc = eccens[i], scale = 0.1)
    return newEccens

planets = pd.read_csv("./data/gasGiantData.csv")

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ") & (planets["pl_sus_data"] == 0)]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_sus_data"] == 0)]
hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0) & (planets["pl_sus_data"] == 0)]

gasGiantEccen = np.array(gasGiants["pl_orbeccen"])
gasGiantEccenErr = (np.array(gasGiants["pl_orbeccenerr1"]) - np.array(gasGiants["pl_orbeccenerr2"]))/2

seCompanionEccen = np.array(seCompanions["pl_orbeccen"])
seCompanionEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

hjCompanionEccen = np.array(hjCompanions["pl_orbeccen"])
hjCompanionEccenErr = (np.array(hjCompanions["pl_orbeccenerr1"]) - np.array(hjCompanions["pl_orbeccenerr2"]))/2

results = np.zeros((5000,3))


for i in range(len(results)):
    gasGiantEccenRedraw = eccenRedraw(gasGiantEccen, gasGiantEccenErr)
    seCompanionEccenRedraw = eccenRedraw(seCompanionEccen, seCompanionEccenErr)
    hjCompanionEccenRedraw = eccenRedraw(hjCompanionEccen, hjCompanionEccenErr)
    ggFitResults = stats.fit(stats.beta, gasGiantEccenRedraw, bounds = ((1,100),(0,100)))
    results[i,0] = ggFitResults.params[1]
    seFitResults = stats.fit(stats.beta, seCompanionEccenRedraw, bounds = ((1,100),(0,100)))
    results[i,1] = seFitResults.params[1]
    hjFitResults = stats.fit(stats.beta, hjCompanionEccenRedraw, bounds = ((1,100),(0,100)))
    results[i,2] = hjFitResults.params[1]


fig = corner.corner(results)
