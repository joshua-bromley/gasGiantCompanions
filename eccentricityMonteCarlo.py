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

gasGiants = planets.loc[((planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")) & (planets["pl_orbeccen"] > 0)]
hotJupiters = planets.loc[(planets["pl_type"] == "HJ") & (planets["pl_orbeccen"] > 0)]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbeccen"] > 0)]
hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0) & (planets["pl_orbeccen"] > 0)]
wcjCompanions = planets.loc[((planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)) & (planets["pl_orbeccen"] > 0)]
loneGasGiants = gasGiants.loc[planets["sy_pnum"] == 1]

gasGiantEccen = np.array(gasGiants["pl_orbeccen"])
gasGiantEccenErr = (np.array(gasGiants["pl_orbeccenerr1"]) - np.array(gasGiants["pl_orbeccenerr2"]))/2

loneGasGiantEccen = np.array(loneGasGiants["pl_orbeccen"])
loneGasGiantEccenErr = (np.array(loneGasGiants["pl_orbeccenerr1"]) - np.array(loneGasGiants["pl_orbeccenerr2"]))/2

hotJupiterEccen = np.array(hotJupiters["pl_orbeccen"])
hotJupiterEccenErr = (np.array(hotJupiters["pl_orbeccenerr1"]) - np.array(hotJupiters["pl_orbeccenerr2"]))/2


seCompanionEccen = np.array(seCompanions["pl_orbeccen"])
seCompanionEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

hjCompanionEccen = np.array(hjCompanions["pl_orbeccen"])
hjCompanionEccenErr = (np.array(hjCompanions["pl_orbeccenerr1"]) - np.array(hjCompanions["pl_orbeccenerr2"]))/2

wcjCompanionEccen = np.array(wcjCompanions["pl_orbeccen"])
wcjCompanionEccenErr = (np.array(wcjCompanions["pl_orbeccenerr1"]) - np.array(wcjCompanions["pl_orbeccenerr2"]))/2

betaResults = np.zeros((2000,8))
rayleighResults = np.zeros((2000,8))

print(len(gasGiantEccen), len(hotJupiterEccen), len(hjCompanionEccen), len(wcjCompanionEccen))

for i in tqdm(range(len(betaResults))):
    gasGiantEccenRedraw = eccenRedraw(loneGasGiantEccen, loneGasGiantEccenErr)
    hotJupiterEccenRedraw = eccenRedraw(hotJupiterEccen, hotJupiterEccenErr)
    hjCompanionEccenRedraw = eccenRedraw(hjCompanionEccen, hjCompanionEccenErr)
    wcjCompanionEccenRedraw = eccenRedraw(wcjCompanionEccen, wcjCompanionEccenErr)
    ggFitResults = stats.fit(stats.beta, gasGiantEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,0] = ggFitResults.params[0]
    betaResults[i,1] = ggFitResults.params[1]
    #rayleighResults[i,0],rayleighResults[i,1] = stats.rayleigh.fit(gasGiantEccenRedraw)
    seFitResults = stats.fit(stats.beta, hotJupiterEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,2] = seFitResults.params[0]
    betaResults[i,3] = seFitResults.params[1]
    #rayleighResults[i,2],rayleighResults[i,3] = stats.rayleigh.fit(seCompanionEccenRedraw)
    hjFitResults = stats.fit(stats.beta, hjCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,4] = hjFitResults.params[0]
    betaResults[i,5] = hjFitResults.params[1]
    rayleighResults[i,4], rayleighResults[i,5] = stats.rayleigh.fit(hjCompanionEccenRedraw)
    wcjFitResults = stats.fit(stats.beta, wcjCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,6] = wcjFitResults.params[0]
    betaResults[i,7] = wcjFitResults.params[1]
    #rayleighResults[i,6],rayleighResults[i,7] = stats.rayleigh.fit(wcjCompanionEccenRedraw)


fig, ax = plt.subplots(1,1)
sb.kdeplot(x = betaResults[:,0], y = betaResults[:,1], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,2], y = betaResults[:,3], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,4], y = betaResults[:,5], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,6], y = betaResults[:,7], ax = ax, levels=[0.01,0.05,0.32])
ax.set_xlabel("$\\alpha$")
ax.set_ylabel("$\\beta$")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:blue", label = "Lone Gas Giants")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:orange", label = "Hot Jupiters")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:green", label = "HJ Companions")
ax.plot(-10,-10, ls = "", marker = ".", color = "tab:red", label = "WCJ Companions")
ax.set_xlim(0.5,3)
ax.set_ylim(1,7)
ax.legend(frameon = False)
fig.savefig("./plots/betaFits3KDE.png")

print(f"Lone Gas Giants a = {np.mean(betaResults[:,0])}+/-{np.std(betaResults[:,0])}, b = {np.mean(betaResults[:,1])}+/-{np.std(betaResults[:,1])}")
print(f"Hot Jupiters a = {np.mean(betaResults[:,2])}+/-{np.std(betaResults[:,2])}, b = {np.mean(betaResults[:,3])}+/-{np.std(betaResults[:,3])}")
print(f"HJ Companions a = {np.mean(betaResults[:,4])}+/-{np.std(betaResults[:,4])}, b = {np.mean(betaResults[:,5])}+/-{np.std(betaResults[:,5])}")
print(f"WCJ Companions a = {np.mean(betaResults[:,6])}+/-{np.std(betaResults[:,6])}, b = {np.mean(betaResults[:,7])}+/-{np.std(betaResults[:,7])}")

print(f"HJ Companions a = {np.mean(rayleighResults[:,4])}+/-{np.std(rayleighResults[:,4])}, b = {np.mean(rayleighResults[:,5])}+/-{np.std(rayleighResults[:,5])}")
'''
fig, ax = plt.subplots(1,1)
sb.kdeplot(x = rayleighResults[:,0], y = rayleighResults[:,1], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = rayleighResults[:,2], y = rayleighResults[:,3], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = rayleighResults[:,4], y = rayleighResults[:,5], ax = ax, levels=[0.01,0.05,0.32])
sb.kdeplot(x = rayleighResults[:,6], y = rayleighResults[:,7], ax = ax, levels=[0.01,0.05,0.32])
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
'''



