import numpy as np
import pandas as pd
from scipy import stats as stats

planets = pd.read_csv("./data/gasGiantData.csv")

gasGiants = planets.loc[((planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")) & (planets["pl_orbeccen"] > 0)]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbeccen"] > 0)]
hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0) & (planets["pl_orbeccen"] > 0)]
wcjCompanions = planets.loc[((planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)) & (planets["pl_orbeccen"] > 0)]

distributions = [gasGiants, seCompanions, hjCompanions, wcjCompanions]
names = ["Gas Giants", "SE Companions", "HJ Companions", "WCJ Companions"]

def chiSq(model, x, y, errs, args):
    chiSq = 0
    for i in range(len(x)):
        chiSq += 0.5*((model(x[i], *args) - y[i])/errs[i])**2
    return chiSq

def eccenRedraw(eccens, errors):
    newEccens = np.zeros_like(eccens)
    for i in range(len(eccens)):
        if not (np.isnan(eccens[i]) or np.isnan(errors[i])):
            while newEccens[i] <= 0:
                newEccens[i] = stats.norm.rvs(loc = eccens[i], scale = errors[i])
        elif np.isnan(errors[i]):
            while newEccens[i] <= 0:
                newEccens[i] = stats.norm.rvs(loc = eccens[i], scale = 0.1)
    return newEccens

def generateCDFerrs(eccens, errors):
    cdfs = np.zeros((10,100))
    for i in range(10):
        if i != 0:
            newEccens = eccenRedraw(eccens, errors)
            newEccens = np.sort(newEccens)
        else:
            newEccens = np.sort(eccens)
        for j in range(100):
            total = 0
            for k in range(len(newEccens)):
                if newEccens[k] < 0.01*j:
                    total += 1
            cdfs[i][j] = total
    
    errors = np.zeros(100)
    for i in range(100):
        errors[i] = np.std(cdfs[:,i])
        if errors[i] == 0:
            errors[i] = 0.01
    return cdfs[0],errors

for i in range(len(distributions)):

    cdf, err = generateCDFerrs(np.array(distributions[i]["pl_orbeccen"]), np.array(distributions[i]["pl_orbeccenerr1"]))
    betaFit = stats.fit(stats.beta, distributions[i]["pl_orbeccen"], bounds = ((0,100),(0,100))).params
    chiSqBeta = chiSq(stats.beta.cdf, np.linspace(0,1,100), cdf, err, betaFit)
    rayleighFit = stats.rayleigh.fit(distributions[i]["pl_orbeccen"])
    chiSqRayleigh = chiSq(stats.rayleigh.cdf, np.linspace(0,1,100), cdf, err, rayleighFit)
    restRayleighFit = stats.rayleigh.fit(distributions[i]["pl_orbeccen"], floc = 0)
    chiSqRestRayleigh = chiSq(stats.rayleigh.cdf, np.linspace(0,1,100), cdf, err, restRayleighFit)
    bicBeta = 2*np.log(len(distributions[i]["pl_orbeccen"])) - 2*np.log(chiSqBeta)
    bicRayleigh = 2*np.log(len(distributions[i]["pl_orbeccen"])) - 2*np.log(chiSqRayleigh)
    bicRestRayleigh = np.log(len(distributions[i]["pl_orbeccen"])) - 2*np.log(chiSqRestRayleigh)

    print(names[i])
    print(f"Beta BIC: {bicBeta}")
    print(f"Rayleigh BIC: {bicRayleigh}")
    print(f"Restricted Rayleigh BIC: {bicRestRayleigh}")


