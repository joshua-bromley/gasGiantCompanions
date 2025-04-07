import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats as stats

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


eccentricityBins = np.arange(0.0,1.01, 0.1)
binCenters = 0.5*(np.add(eccentricityBins[:-1], eccentricityBins[1:]))

eccenDist = np.zeros((len(gasGiantsList), len(binCenters), 3))

for i in range(len(gasGiantsList)):
    gasGiants = gasGiantsList[i]
    nPlanets = len(gasGiants)
    planetCounts = np.zeros(len(binCenters))
    for j in range(len(gasGiants)):
        eccen = gasGiants.iloc[j]["pl_orbeccen"]
        binIdx = int(eccen*10)
        planetCounts[binIdx] += 1
    
    for j in range(len(planetCounts)):
        a = planetCounts[j] + 1
        b = nPlanets - planetCounts[j] + 1
        eccenDist[i][j][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.68, a, b)
        eccenDist[i][j][1] = eccenDist[i][j][0] - errorBars[0]
        eccenDist[i][j][2] = errorBars[1] - eccenDist[i][j][0]


fig, ax = plt.subplots(1,1)

for i in range(len(eccenDist)):
    ax.errorbar(binCenters, eccenDist[i][:,0], yerr = eccenDist[i][:,1:].T, ls = "", marker = symbols[i], label = gasGiantLabels[i], alpha = 0.8, capsize = 5)

ax.legend(frameon = False)
ax.set_ylabel("Occurrence")
ax.set_xlabel("Eccentricity")
plt.show()

