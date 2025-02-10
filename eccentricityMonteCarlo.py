import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats
import corner
from tqdm import tqdm
import seaborn as sb

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

planets = pd.read_csv("./data/gasGiantData.csv")

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ") & (planets["pl_orbeccen"] > 0)]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbeccen"] > 0)]
hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0) & (planets["pl_orbeccen"] > 0)]
wcjCompanions = planets.loc[((planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)) & (planets["pl_orbeccen"] > 0)]

gasGiantEccen = np.array(gasGiants["pl_orbeccen"])
gasGiantEccenErr = (np.array(gasGiants["pl_orbeccenerr1"]) - np.array(gasGiants["pl_orbeccenerr2"]))/2

seCompanionEccen = np.array(seCompanions["pl_orbeccen"])
seCompanionEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

hjCompanionEccen = np.array(hjCompanions["pl_orbeccen"])
hjCompanionEccenErr = (np.array(hjCompanions["pl_orbeccenerr1"]) - np.array(hjCompanions["pl_orbeccenerr2"]))/2

wcjCompanionEccen = np.array(wcjCompanions["pl_orbeccen"])
wcjCompanionEccenErr = (np.array(wcjCompanions["pl_orbeccenerr1"]) - np.array(wcjCompanions["pl_orbeccenerr2"]))/2

betaResults = np.zeros((2000,8))
rayleighResults = np.zeros((2000,8))


for i in tqdm(range(len(betaResults))):
    gasGiantEccenRedraw = eccenRedraw(gasGiantEccen, gasGiantEccenErr)
    seCompanionEccenRedraw = eccenRedraw(seCompanionEccen, seCompanionEccenErr)
    hjCompanionEccenRedraw = eccenRedraw(hjCompanionEccen, hjCompanionEccenErr)
    wcjCompanionEccenRedraw = eccenRedraw(wcjCompanionEccen, wcjCompanionEccenErr)
    ggFitResults = stats.fit(stats.beta, gasGiantEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,0] = ggFitResults.params[0]
    betaResults[i,1] = ggFitResults.params[1]
    rayleighResults[i,0],rayleighResults[i,1] = stats.rayleigh.fit(gasGiantEccenRedraw)
    seFitResults = stats.fit(stats.beta, seCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,2] = seFitResults.params[0]
    betaResults[i,3] = seFitResults.params[1]
    rayleighResults[i,2],rayleighResults[i,3] = stats.rayleigh.fit(seCompanionEccenRedraw)
    hjFitResults = stats.fit(stats.beta, hjCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,4] = hjFitResults.params[0]
    betaResults[i,5] = hjFitResults.params[1]
    rayleighResults[i,4], rayleighResults[i,5] = stats.rayleigh.fit(hjCompanionEccenRedraw)
    wcjFitResults = stats.fit(stats.beta, wcjCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,6] = wcjFitResults.params[0]
    betaResults[i,7] = wcjFitResults.params[1]
    rayleighResults[i,6],rayleighResults[i,7] = stats.rayleigh.fit(wcjCompanionEccenRedraw)


fig, ax = plt.subplots(1,1)
sb.kdeplot(x = betaResults[:,0], y = betaResults[:,1], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = betaResults[:,2], y = betaResults[:,3], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = betaResults[:,4], y = betaResults[:,5], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = betaResults[:,6], y = betaResults[:,7], ax = ax, levels=[0.68,0.95,0.99])
ax.set_xlabel("$\\alpha$")
ax.set_ylabel("$\\beta$")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:blue", label = "All Gas Giants")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:orange", label = "SE Companions")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:green", label = "HJ Companions")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:red", label = "WCJ Companions")
ax.set_xlim(0.5,3)
ax.set_ylim(1,7)
ax.legend(frameon = False)
fig.savefig("./plots/betaFitsKDE.png")

fig, ax = plt.subplots(1,1)
sb.kdeplot(x = rayleighResults[:,0], y = rayleighResults[:,1], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = rayleighResults[:,2], y = rayleighResults[:,3], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = rayleighResults[:,4], y = rayleighResults[:,5], ax = ax, levels=[0.68,0.95,0.99])
sb.kdeplot(x = rayleighResults[:,6], y = rayleighResults[:,7], ax = ax, levels=[0.68,0.95,0.99])
ax.set_xlabel("loc")
ax.set_ylabel("$\sigma$")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:blue", label = "All Gas Giants")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:orange", label = "SE Companions")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:green", label = "HJ Companions")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:red", label = "WCJ Companions")
ax.legend(frameon = False)
ax.set_xlim(-0.15,0.1)
ax.set_ylim(0.15,0.4)
fig.savefig("./plots/rayleighFitsKDE.png")



