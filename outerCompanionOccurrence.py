import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats
import pickle as pk

clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows=101)
clsStars = pd.read_csv("./data/clsStars.csv")
planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[planets["baseline"] > 1*365]

clsGasGiants = clsPlanets.loc[(clsPlanets["pl_bmassj"] > 0.5) & (clsPlanets["pl_orbsmax"] > 0.1)]
#clsStars = clsStars.loc[clsStars["[Fe/H]"] > 0.041]
#planets = planets.loc[planets["st_met"] > 0.041]

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
hotJupiters = planets.loc[planets["pl_type"] == "HJ"]
subSaturns = planets.loc[(planets["pl_type"] == "HS") | (planets["pl_type"] == "CS")]
superEarths = planets.loc[planets["pl_type"] == "SE"]

seCompanions = planets.loc[planets["companion_type"] % 2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
wcjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]

def occurrence(x, a, b):
    return (1/sp.beta(a,b)) * x**(a-1) * (1-x)**(b-1)

planetPopulations = [clsGasGiants, hjCompanions, wcjCompanions]
hostPopulations = [clsStars, hotJupiters, gasGiants]
labels = ["P(GG)", "P(GG|HJ)", "P(GG|WCJ)"]

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
    
    a = n_det -1
    b = n_tot- n_missed - n_det -1 
    
    occurrenceRates[i][0] = stats.beta.median(a,b)
    errorBars = stats.beta.interval(0.68, a, b)
    occurrenceRates[i][1] = occurrenceRates[i][0] - errorBars[0]
    occurrenceRates[i][2] = errorBars[1] - occurrenceRates[i][0]
    ax.plot(x, occurrence(x,a,b), label = labels[i])
    if i > 0:
        enhancement[i] = (occurrenceRates[i][0] - occurrenceRates[0][0])/(np.sqrt(occurrenceRates[i][1]**2 + occurrenceRates[0][2]**2))

    print(labels[i])
    print(f"{n_det}/{n_tot-n_missed}, {occurrenceRates[i][0]} (-{occurrenceRates[i][1]} +{occurrenceRates[i][2]}), {enhancement[i]}")
    print(n_tot, n_missed)

ax.legend(frameon = False)
ax.set_xlabel("P(GG)")
ax.set_ylabel("PDF")
#fig.savefig("./plots/outerCompanionOccurrence.png")

plt.show()


