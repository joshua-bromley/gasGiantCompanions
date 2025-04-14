import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib  import rcParams
from scipy import special as sp
from scipy import stats as stats
import corner
from tqdm import tqdm
import seaborn as sb

##Adjust plotting defaults
rcParams["axes.linewidth"] = 2

rcParams["ytick.right"] = True
rcParams["ytick.direction"] = "in"
rcParams["ytick.minor.visible"] = True
rcParams["ytick.major.left"] = True
rcParams["ytick.major.right"] = True
rcParams["ytick.minor.left"] = True
rcParams["ytick.minor.right"] = True
rcParams["ytick.major.size"] = 10
rcParams["ytick.minor.size"] = 5
rcParams["ytick.major.width"] = 1
rcParams["ytick.minor.width"] = 1


rcParams["xtick.top"] = True
rcParams["xtick.direction"] = "in"
rcParams["xtick.minor.visible"] = True
rcParams["xtick.major.top"] = True
rcParams["xtick.major.bottom"] = True
rcParams["xtick.minor.top"] = True
rcParams["xtick.minor.bottom"] = True
rcParams["xtick.major.size"] = 10
rcParams["xtick.minor.size"] = 5
rcParams["xtick.major.width"] = 1
rcParams["xtick.minor.width"] = 1

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

planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]


gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
cjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]

gasGiantEccen = np.array(gasGiants["pl_orbeccen"])
gasGiantEccenErr = (np.array(gasGiants["pl_orbeccenerr1"]) - np.array(gasGiants["pl_orbeccenerr2"]))/2

#loneGasGiantEccen = np.array(loneGasGiants["pl_orbeccen"])
#loneGasGiantEccenErr = (np.array(loneGasGiants["pl_orbeccenerr1"]) - np.array(loneGasGiants["pl_orbeccenerr2"]))/2

#hotJupiterEccen = np.array(hotJupiters["pl_orbeccen"])
#hotJupiterEccenErr = (np.array(hotJupiters["pl_orbeccenerr1"]) - np.array(hotJupiters["pl_orbeccenerr2"]))/2


seCompanionEccen = np.array(seCompanions["pl_orbeccen"])
seCompanionEccenErr = (np.array(seCompanions["pl_orbeccenerr1"]) - np.array(seCompanions["pl_orbeccenerr2"]))/2

hjCompanionEccen = np.array(hjCompanions["pl_orbeccen"])
hjCompanionEccenErr = (np.array(hjCompanions["pl_orbeccenerr1"]) - np.array(hjCompanions["pl_orbeccenerr2"]))/2

cjCompanionEccen = np.array(cjCompanions["pl_orbeccen"])
cjCompanionEccenErr = (np.array(cjCompanions["pl_orbeccenerr1"]) - np.array(cjCompanions["pl_orbeccenerr2"]))/2

ssCompanionEccen = np.array(ssCompanions["pl_orbeccen"])
ssCompanionEccenErr = (np.array(ssCompanions["pl_orbeccenerr1"]) - np.array(ssCompanions["pl_orbeccenerr2"]))/2


betaResults = np.zeros((2000,10))
rayleighResults = np.zeros((2000,10))

print(len(gasGiantEccen), len(seCompanionEccen), len(hjCompanionEccen), len(cjCompanionEccen))

for i in tqdm(range(len(betaResults))):
    gasGiantEccenRedraw = eccenRedraw(gasGiantEccen, gasGiantEccenErr)
    seCompanionEccenRedraw = eccenRedraw(seCompanionEccen, seCompanionEccenErr)
    hjCompanionEccenRedraw = eccenRedraw(hjCompanionEccen, hjCompanionEccenErr)
    wcjCompanionEccenRedraw = eccenRedraw(cjCompanionEccen, cjCompanionEccenErr)
    ssCompanionEccenRedraw = eccenRedraw(ssCompanionEccen, ssCompanionEccenErr)
    ggFitResults = stats.fit(stats.beta, gasGiantEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,0] = ggFitResults.params[0]
    betaResults[i,1] = ggFitResults.params[1]
    #rayleighResults[i,0],rayleighResults[i,1] = stats.rayleigh.fit(gasGiantEccenRedraw)
    seFitResults = stats.fit(stats.beta, seCompanionEccenRedraw, bounds = ((0,100),(0,100)))
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
    ssFitResults = stats.fit(stats.beta, ssCompanionEccenRedraw, bounds = ((0,100),(0,100)))
    betaResults[i,8] = ssFitResults.params[0]
    betaResults[i,9] = ssFitResults.params[1]
    #rayleighResults[i,6],rayleighResults[i,7] = stats.rayleigh.fit(wcjCompanionEccenRedraw)


fig, ax = plt.subplots(1,2, figsize = (10,4))
sb.kdeplot(x = betaResults[:,0], y = betaResults[:,1], ax = ax[0], levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,2], y = betaResults[:,3], ax = ax[0], levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,8], y = betaResults[:,9], ax = ax[0], levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,4], y = betaResults[:,5], ax = ax[0], levels=[0.01,0.05,0.32])
sb.kdeplot(x = betaResults[:,6], y = betaResults[:,7], ax = ax[0], levels=[0.01,0.05,0.32])

ax[0].set_xlabel("$\\alpha$", fontsize = 16)
ax[0].set_ylabel("$\\beta$", fontsize = 16)
ax[0].plot(-10,-10, ls = "", marker = ".", color = "tab:blue", label = "Cold Jupiters")
ax[0].plot(-10,-10, ls = "", marker = ".", color = "tab:orange", label = "SE Companions")
ax[0].plot(-10,-10, ls = "", marker = ".", color = "tab:green", label = "SS Companions")
ax[0].plot(-10,-10, ls = "", marker = ".", color = "tab:red", label = "HJ Companions")
ax[0].plot(-10,-10, ls = "", marker = ".", color = "tab:purple", label = "CJ Companions")

xTicks = np.arange(0.5,3.1,0.5)
yTicks = np.arange(1,7.1,1)
ax[0].set_xticks(xTicks)
ax[0].set_yticks(yTicks)
ax[0].set_xlim(0.5,3)
ax[0].set_ylim(1,7)
ax[0].legend(frameon = False)

tickLabelSize = 12
ax[0].tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)


x = np.linspace(0,1,100)
ax[1].plot(x, stats.beta.pdf(x,np.mean(betaResults[:,0]),np.mean(betaResults[:,1])), color = "tab:blue")
ax[1].plot(x, stats.beta.pdf(x,np.mean(betaResults[:,2]),np.mean(betaResults[:,3])), color = "tab:orange")
ax[1].plot(x, stats.beta.pdf(x,np.mean(betaResults[:,4]),np.mean(betaResults[:,5])), color = "tab:red")
ax[1].plot(x, stats.beta.pdf(x,np.mean(betaResults[:,6]),np.mean(betaResults[:,7])), color = "tab:purple")
ax[1].plot(x, stats.beta.pdf(x,np.mean(betaResults[:,8]),np.mean(betaResults[:,9])), color = "tab:green")

ax[1].set_xlabel("Eccentricity", fontsize = 16)
ax[1].set_ylabel("Probability Density", fontsize = 16)

xTicks = np.arange(0,1.01,0.2)
yTicks = np.arange(0,4.1,0.5)
ax[1].set_xticks(xTicks)
ax[1].set_yticks(yTicks)
ax[1].set_xlim(-0.05,1.05)
ax[1].set_ylim(0.0,4.0)

ax[1].tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()


fig.savefig("./plots/betaFits4KDE.png")
fig.savefig("./plots/betaFits4KDE.pdf")

print(f"Gas Giants a = {np.mean(betaResults[:,0])}+/-{np.std(betaResults[:,0])}, b = {np.mean(betaResults[:,1])}+/-{np.std(betaResults[:,1])}")
print(f"SE Companions a = {np.mean(betaResults[:,2])}+/-{np.std(betaResults[:,2])}, b = {np.mean(betaResults[:,3])}+/-{np.std(betaResults[:,3])}")
print(f"HJ Companions a = {np.mean(betaResults[:,4])}+/-{np.std(betaResults[:,4])}, b = {np.mean(betaResults[:,5])}+/-{np.std(betaResults[:,5])}")
print(f"WCJ Companions a = {np.mean(betaResults[:,6])}+/-{np.std(betaResults[:,6])}, b = {np.mean(betaResults[:,7])}+/-{np.std(betaResults[:,7])}")
print(f"SS Companions a = {np.mean(betaResults[:,8])}+/-{np.std(betaResults[:,8])}, b = {np.mean(betaResults[:,9])}+/-{np.std(betaResults[:,9])}")


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



