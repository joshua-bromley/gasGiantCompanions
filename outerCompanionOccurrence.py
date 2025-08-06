import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats
import pickle as pk
from matplotlib import rcParams

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

clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows=101)
clsStars = pd.read_csv("./data/clsStars.csv")
planets = pd.read_csv("./data/gasGiantComplete2.csv")
#planets = planets.loc[planets["baseline"] > 1*365]

#clsPlanets = clsPlanets.loc[clsPlanets["st_met"] < 0]
#clsStars = clsStars.loc[clsStars["[Fe/H]"] < 0]
#planets = planets.loc[planets["st_met"] < 0]

clsGasGiants = clsPlanets.loc[(clsPlanets["pl_bmassj"] > 0.5) & (clsPlanets["pl_orbsmax"] > 0.1)]
#clsStars = clsStars.loc[clsStars["[Fe/H]"] > 0.041]
#planets = planets.loc[planets["st_met"] > 0.041]

gasGiants = planets.loc[(planets["pl_type"] == "CJ")]
warmJupiters = planets.loc[(planets["pl_type"] == "WJ") ]
hotJupiters = planets.loc[planets["pl_type"] == "HJ"]
hotSaturns = planets.loc[(planets["pl_type"] == "HS")]
subSaturns = planets.loc[(planets["pl_type"] == "CS")]
superEarths = planets.loc[planets["pl_type"] == "SE"]

seCompanions = planets.loc[planets["companion_type"] % 2 == 0]
hsCompanions = planets.loc[(planets["companion_type"] % 3 == 0) ]
csCompanions = planets.loc[(planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
wjCompanions = planets.loc[(planets["companion_type"] % 11 == 0)]
cjCompanions = planets.loc[(planets["companion_type"] % 13 == 0)]

def occurrence(x, a, b):
    return (1/sp.beta(a,b)) * x**(a-1) * (1-x)**(b-1)

planetPopulations = [clsGasGiants, seCompanions, hsCompanions, csCompanions, hjCompanions, wjCompanions, cjCompanions]
hostPopulations = [clsStars, superEarths, hotSaturns, subSaturns, hotJupiters, warmJupiters, gasGiants]
labels = ["P(GG) (R21)", "P(GG|SE)", "P(GG|HS)", "P(GG|CS)", "P(GG|HJ)", "P(GG|WJ)", "P(GG|CJ)"]
colours = ["tab:blue", "tab:orange", "tab:green", "tab:cyan", "tab:red", "tab:pink", "tab:purple"]

occurrenceRates = np.zeros((len(planetPopulations),3))
enhancement = np.zeros(len(planetPopulations))

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

fig, ax = plt.subplots(1,1, figsize = (6,4))

x = np.linspace(0,0.5,1000)


for i in range(len(planetPopulations)):
    n_det = len(np.unique(planetPopulations[i]["hostname"]))
    n_missed = 0
    if i > 0:
        n_tot = len(np.unique(hostPopulations[i]["hostname"]))
        hostPopulation = hostPopulations[i].drop_duplicates(subset="hostname")
        for j in range(len(hostPopulation)):
            hostname = hostPopulation.iloc[j]["hostname"]
            compProb = pk.load(open("./completenessMaps/"+hostname+".p", "rb"))
            probs = compProb["prob"]
            tot = ((maxSmaxIdx - minSmaxIdx)*(maxMassIdx - minMassIdx))#(30.*40.)#
            sumProb = np.sum(probs[minSmaxIdx:maxSmaxIdx,minMassIdx:maxMassIdx])#[0:30,5:45])#
            n_missed += (tot-sumProb)/tot
            #if i == 1:
                #print(hostname, (tot-sumProb)/tot)
    else:
        n_tot = len(np.unique(hostPopulations[i]["CPS identifier"]))
    
    a = n_det +1
    b = n_tot- n_missed - n_det +1 
    
    occurrenceRates[i][0] = stats.beta.median(a,b)
    errorBars = stats.beta.interval(0.68, a, b)
    occurrenceRates[i][1] = occurrenceRates[i][0] - errorBars[0]
    occurrenceRates[i][2] = errorBars[1] - occurrenceRates[i][0]
    ax.plot(x, occurrence(x,a,b), label = labels[i], color = colours[i])
    if i > 0:
        enhancement[i] = (occurrenceRates[i][0] - occurrenceRates[0][0])/(np.sqrt(occurrenceRates[i][1]**2 + occurrenceRates[0][2]**2))

    print(labels[i])
    print(f"{n_det}/{n_tot-n_missed}, {occurrenceRates[i][0]} (-{occurrenceRates[i][1]} +{occurrenceRates[i][2]}), {enhancement[i]}")
    print(n_tot, n_missed)

ax.legend(frameon = False, fontsize = 12)
ax.set_xlabel("P(GG)", fontsize = 16)
ax.set_ylabel("PDF", fontsize = 16)

yTicks = np.arange(0,30.1,5)
xTicks = np.arange(0,0.51,0.1)

ax.set_xticks(xTicks)
ax.set_yticks(yTicks)
ax.set_xlim(-0.02,0.52)
ax.set_ylim(-1,31)

#ax.text(0.2,27.5, "$[Fe/H] < 0$ (Metal Poor)", fontsize = 18)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()

#fig.savefig("./plots/outerCompanionOccurrenceMetalPoor.png")

plt.show()


