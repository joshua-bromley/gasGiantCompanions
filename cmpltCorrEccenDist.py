import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats as stats
import pickle as pk
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

planets = pd.read_csv("./data/gasGiantComplete2.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]


coldJupiters = planets.loc[ (planets["pl_type"] == "CJ")]
warmJupiters = planets.loc[(planets["pl_type"] == "WJ")]
hotJupiters = planets.loc[planets["pl_type"] == "HJ"]
superEarths = planets.loc[planets["pl_type"] == "SE"]
hotSaturns = planets.loc[planets["pl_type"] == "HS"]
coldSaturns = planets.loc[planets["pl_type"] == "CS"]

seCompanions = planets.loc[planets["companion_type"] %2 == 0]
hsCompanions = planets.loc[(planets["companion_type"] % 3 == 0)]
csCompanions = planets.loc[(planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
wjCompanions = planets.loc[(planets["companion_type"] % 11 == 0)]
cjCompanions = planets.loc[(planets["companion_type"] % 13 == 0)]

gasGiantsList = [ hjCompanions, cjCompanions, seCompanions, wjCompanions, hsCompanions, csCompanions]
hostPopulations = [hotJupiters, coldJupiters, superEarths, warmJupiters, hotSaturns, coldSaturns]
gasGiantLabels= [ "HJ Companions", "CJ Companions", "SE Companions", "WJ Companions", "HS Companions", "CS Companions"]
symbols = [ '*', "P", "s", "x", "^", "D"]
colours = [ "tab:red", "tab:purple", "tab:orange", "tab:pink", "tab:green", "tab:cyan"]


cmpltMapMassBins = np.logspace(np.log10(0.3), np.log10(30), 50)
cmpltMapSmaxBins = np.logspace(np.log10(0.3), np.log10(30), 50)
minMass = 0.5
maxMass = 20
minSmax = 1
maxSmax = 10
minSmaxIdx = int(np.abs(minSmax - cmpltMapSmaxBins).argmin())
maxSmaxIdx = int(np.abs(maxSmax - cmpltMapSmaxBins).argmin())
minMassIdx = int(np.abs(minMass - cmpltMapMassBins).argmin())
maxMassIdx = int(np.abs(maxMass - cmpltMapMassBins).argmin())


betaDistParams = [[1.07,0.78,1.18,1.90,1.22],[3.16,2.03,3.97,2.76,4.55]]
rayleighModes = [0.265,0.306,0.244,0.330,0.229]


eccentricityBins = np.arange(0.0,1.01, 0.2)
binCenters = 0.5*(np.add(eccentricityBins[:-1], eccentricityBins[1:]))
eccenBinWidth = np.mean(np.subtract(eccentricityBins[1:],eccentricityBins[:-1] ))

eccenDist = np.zeros((len(gasGiantsList), len(binCenters), 3))

for i in range(len(gasGiantsList)):
    gasGiants = gasGiantsList[i]
    nPlanets = len(hostPopulations[i])
    planetCounts = np.zeros(len(binCenters))
    for j in range(len(gasGiants)):
        eccen = gasGiants.iloc[j]["pl_orbeccen"]
        binIdx = int(eccen*len(binCenters))
        planetCounts[binIdx] += 1
    if i == 3:
        print(planetCounts)
    n_missed = 0
    for j in range(len(hostPopulations[i])):
            hostname = hostPopulations[i].iloc[j]["hostname"]
            compProb = pk.load(open("./completenessMaps/"+hostname+".p", "rb"))
            probs = compProb["prob"]
            tot = ((maxSmaxIdx - minSmaxIdx)*(maxMassIdx - minMassIdx))#(30.*40.)#
            sumProb = np.sum(probs[minSmaxIdx:maxSmaxIdx,minMassIdx:maxMassIdx])#[0:30,5:45])#
            n_missed += (tot-sumProb)/tot
    
    for j in range(len(planetCounts)):
        a = planetCounts[j] + 1
        b = nPlanets - n_missed - planetCounts[j] + 1
        eccenDist[i][j][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.95, a, b)
        eccenDist[i][j][1] = eccenDist[i][j][0] - errorBars[0]
        eccenDist[i][j][2] = errorBars[1] - eccenDist[i][j][0]


fig, ax = plt.subplots(1,1)
x = np.linspace(0,1,1000)

for i in (5,2,1,0):
    ax.errorbar(binCenters, eccenDist[i][:,0]/eccenBinWidth, yerr = eccenDist[i][:,1:].T/eccenBinWidth, ls = "", marker = symbols[i], label = gasGiantLabels[i], alpha = 0.8, capsize = 5, color = colours[i])
    #ax.plot(x, stats.beta.pdf(x, a = betaDistParams[0][i], b = betaDistParams[1][i]), color = colours[i], alpha = 0.5)
    #ax.plot(x, stats.rayleigh.pdf(x, scale = rayleighModes[i]), color = colours[i], alpha = 0.5)


xTicks = np.arange(0.0,1.01,0.2)
yTicks = np.arange(0.0,1.1,0.2)

ax.set_xticks(xTicks)
ax.set_yticks(yTicks)
ax.set_xlim(0,1)
ax.set_ylim(-0.05,1.05)

ax.legend(frameon = False, fontsize = 12)
ax.set_ylabel("Occurrence", fontsize = 16)
ax.set_xlabel("Eccentricity", fontsize = 16)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()
fig.savefig("./plots/compltCorrEccenDist.png")
#ig.savefig("./plots/compltCorrEccenDist.pdf")
plt.show()

