import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats as stats
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

planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]


coldJupiters = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
cjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]

gasGiantsList = [coldJupiters, seCompanions, ssCompanions, hjCompanions, cjCompanions]
gasGiantLabels= ["Cold Jupiters", 'SE Companions', "SS Companions",  "HJ Companions", "CJ Companions"]
symbols = ["o", "s", 'h', '*', "P"]
colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

betaDistParams = [[1.07,0.78,1.18,1.90,1.22],[3.16,2.03,3.97,2.76,4.55]]
rayleighModes = [0.265,0.306,0.244,0.330,0.229]


eccentricityBins = np.arange(0.0,1.01, 0.2)
binCenters = 0.5*(np.add(eccentricityBins[:-1], eccentricityBins[1:]))
eccenBinWidth = np.mean(np.subtract(eccentricityBins[1:],eccentricityBins[:-1] ))

eccenDist = np.zeros((len(gasGiantsList), len(binCenters), 3))

for i in range(len(gasGiantsList)):
    gasGiants = gasGiantsList[i]
    nPlanets = len(gasGiants)
    planetCounts = np.zeros(len(binCenters))
    for j in range(len(gasGiants)):
        eccen = gasGiants.iloc[j]["pl_orbeccen"]
        binIdx = int(eccen*len(binCenters))
        planetCounts[binIdx] += 1
    if i == 3:
        print(planetCounts)
    
    for j in range(len(planetCounts)):
        a = planetCounts[j] + 1
        b = nPlanets - planetCounts[j] + 1
        eccenDist[i][j][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.95, a, b)
        eccenDist[i][j][1] = eccenDist[i][j][0] - errorBars[0]
        eccenDist[i][j][2] = errorBars[1] - eccenDist[i][j][0]


fig, ax = plt.subplots(1,1)
x = np.linspace(0,1,1000)

for i in range(len(eccenDist)):
    ax.errorbar(binCenters, eccenDist[i][:,0]/eccenBinWidth, yerr = eccenDist[i][:,1:].T/eccenBinWidth, ls = "", marker = symbols[i], label = gasGiantLabels[i], alpha = 0.8, capsize = 5)
    ax.plot(x, stats.beta.pdf(x, a = betaDistParams[0][i], b = betaDistParams[1][i]), color = colours[i], alpha = 0.5)
    #ax.plot(x, stats.rayleigh.pdf(x, scale = rayleighModes[i]), color = colours[i], alpha = 0.5)


xTicks = np.arange(0.0,1.01,0.2)
yTicks = np.arange(0.0,4.1,1)

ax.set_xticks(xTicks)
ax.set_yticks(yTicks)
ax.set_xlim(0,1)
ax.set_ylim(-0.2,4.2)

ax.legend(frameon = False)
ax.set_ylabel("Occurrence", fontsize = 16)
ax.set_xlabel("Eccentricity", fontsize = 16)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()
#fig.savefig("./plots/compltCorrEccenDist.png")
#ig.savefig("./plots/compltCorrEccenDist.pdf")
plt.show()

