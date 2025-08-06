import numpy as np
import pandas as pd
from scipy import stats as stats
import pickle as pk

planets = pd.read_csv("./data/gasGiantComplete2.csv")

smallPlanets = planets.loc[planets["pl_type"] == "SE"]
superEarths = smallPlanets.loc[(smallPlanets["pl_bmasse"] < 5)]
miniNeptunes = smallPlanets.loc[(smallPlanets["pl_bmasse"] >= 5)]

smallPlanetCompanions = planets.loc[planets["companion_type"] % 2 == 0]

companionType = np.ones(len(smallPlanetCompanions))

for i in range(len(smallPlanetCompanions)):
    hostname = smallPlanetCompanions.iloc[i]["hostname"]
    seInner = superEarths.loc[superEarths["hostname"] == hostname]
    mnINner = miniNeptunes.loc[miniNeptunes["hostname"] == hostname]
    if len(seInner) > 0:
        companionType[i] *= 2
    if len(mnINner) > 0:
        companionType[i] *= 3

smallPlanetCompanions.insert(20, "seCompanionType", companionType)

mnCompanions = smallPlanetCompanions.loc[smallPlanetCompanions["seCompanionType"] % 3 == 0]
seCompanions = smallPlanetCompanions.loc[smallPlanetCompanions["seCompanionType"] % 2 == 0]

print(len(seCompanions), len(superEarths))
print(len(mnCompanions), len(miniNeptunes))
print(len(smallPlanetCompanions), len(smallPlanets))



planetPopulations = [mnCompanions, seCompanions]
hostPopulations = [miniNeptunes, superEarths]
labels = ["P(GG|MN)", "P(GG|SE)"]

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

for i in range(len(planetPopulations)):
    n_det = len(np.unique(planetPopulations[i]["hostname"]))
    n_missed = 0
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

    
    a = n_det +1
    b = n_tot- n_missed - n_det +1 
    
    occurrenceRates[i][0] = stats.beta.median(a,b)
    errorBars = stats.beta.interval(0.68, a, b)
    occurrenceRates[i][1] = occurrenceRates[i][0] - errorBars[0]
    occurrenceRates[i][2] = errorBars[1] - occurrenceRates[i][0]
    #ax.plot(x, occurrence(x,a,b), label = labels[i], color = colours[i])
    if i > 0:
        enhancement[i] = (occurrenceRates[i][0] - occurrenceRates[0][0])/(np.sqrt(occurrenceRates[i][1]**2 + occurrenceRates[0][2]**2))

    print(labels[i])
    print(f"{n_det}/{n_tot-n_missed}, {occurrenceRates[i][0]} (-{occurrenceRates[i][1]} +{occurrenceRates[i][2]}), {enhancement[i]}")
    print(n_tot, n_missed)