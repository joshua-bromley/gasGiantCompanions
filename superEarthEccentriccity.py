import numpy as np
import pandas as pd
from scipy import stats as stats
from scipy import special as sp
from matplotlib import pyplot as plt

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

def occurrence(x, a, b):
    return (1/sp.beta(a,b))*x**(a-1) * (1-x)**(b-1)


clsStars = pd.read_csv("./data/clsStars.csv")
clsStars = clsStars.loc[(clsStars["mass"] >= 0.6) & (clsStars["[Fe/H]"] > 0)]
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmasse"] > 0.5*317) & (clsPlanets["pl_orbsmax"] > 1) & (clsPlanets["st_mass"] >= 0.6) & (clsPlanets["st_met"] > 0)]

planets = pd.read_csv("./data/gasGiantData.csv")
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["st_met"] > 0)]
superEarths = planets.loc[(planets["pl_type"] == "SE") & (planets["st_met"] > 0)]

seCompanionsEccen = seCompanions["pl_orbeccen"].dropna().values
seCompanionsEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

seCompanionsRedraw = eccenRedraw(seCompanionsEccen, seCompanionsEccenErr)

fig, ax = plt.subplots(1,1) 
ax.ecdf(seCompanionsEccen, label  = "SE Companions")
ax.ecdf(seCompanionsRedraw)
ax.ecdf(clsPlanets["pl_orbeccen"].dropna(), label = "CLS")
ax.legend(frameon = False)
#plt.show()

nDynHotSECompanions = len(np.unique(seCompanions.loc[seCompanions["pl_orbeccen"] > 0.3]["hostname"]))
nDynColdSECompanions = len(np.unique(seCompanions.loc[seCompanions["pl_orbeccen"] <= 0.3]["hostname"]))
nDynHotCLS = np.count_nonzero(clsPlanets["pl_orbeccen"].values > 0.3)
nDynColdClS = np.count_nonzero(clsPlanets["pl_orbeccen"].values <= 0.3)

nCLSStars = len(clsStars)
nSESystems = len(np.unique(superEarths["hostname"]))



redrawDict = {"hostname": seCompanions.loc[seCompanions["pl_orbeccen"] > -0.00001]["hostname"].values, "pl_orbeccen":seCompanionsRedraw}
seCompanionsRD = pd.DataFrame(redrawDict)

nDynHotSECompanionsRD = len(np.unique(seCompanionsRD.loc[seCompanionsRD["pl_orbeccen"] > 0.3]["hostname"]))
nDynColdSECompanionsRD = len(np.unique(seCompanionsRD.loc[seCompanionsRD["pl_orbeccen"] <= 0.3]["hostname"]))

#a = n_det + 1
#b = n_tot - n_det + 1
x = np.linspace(0,0.5,1000)
fig, ax = plt.subplots(1,1)
ax.plot(x, occurrence(x, nDynHotSECompanions + 1, nSESystems - nDynHotSECompanions + 1), color = "tab:red", ls = "--", label = "P(GG, e>0.3|SE, [Fe/H] > 0)")
ax.plot(x, occurrence(x, nDynHotCLS + 1, nCLSStars - nDynHotCLS + 1), color = "tab:red", label = "P(GG,e>0.3|[Fe/H] > 0)")
ax.plot(x, occurrence(x, nDynHotSECompanionsRD + 1, nSESystems - nDynHotSECompanionsRD + 1), color = "tab:red", ls = "-.")
ax.plot(x, occurrence(x, nDynColdSECompanions + 1, nSESystems - nDynColdSECompanions + 1), color = "tab:blue", ls = "--", label = "P(GG,e<=0.3|SE,[Fe/H]>0)")
ax.plot(x, occurrence(x, nDynColdClS + 1, nCLSStars - nDynColdClS + 1), color = "tab:blue", label = "P(GG,e<=0.3|[Fe/H] > 0)")
ax.plot(x, occurrence(x, nDynColdSECompanionsRD + 1, nSESystems - nDynColdSECompanionsRD + 1), color = "tab:blue", ls = "-.")
ax.legend(frameon = False)
plt.show()