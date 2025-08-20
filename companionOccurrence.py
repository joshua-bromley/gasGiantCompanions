import numpy as np
import pandas as pd
import pickle as pk
from scipy import stats as stats
from matplotlib import pyplot as plt

def occurrenceRate(targetPopulation, hostPopulation, minMassIdx, maxMassIdx, minSmaxIdx, maxSmaxIdx, clsFlag = False):
    n_det = len(np.unique(targetPopulation["hostname"]))
    n_missed = 0

    n_tot = len(np.unique(hostPopulation["hostname"]))
    hostPopulation = hostPopulation.drop_duplicates(subset="hostname")
    if not clsFlag:
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
        n_missed = 0

    
    a = n_det +1
    b = n_tot- n_missed - n_det +1 
    if n_det > 0:
        occurRate = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.68, a, b)
    else:
        errorBars = (0, stats.beta.median(a,b))
        occurRate = 0#(errorBars[0] + errorBars[1])/2
    print(n_det, n_tot, occurRate)
    return np.array((occurRate, occurRate - errorBars[0], errorBars[1] - occurRate)), np.array((n_det, n_tot, n_missed))

planets = pd.read_csv("./data/gasGiantDataComplete2.csv")
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows=101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmassj"] > 0.5) & (clsPlanets["pl_orbsmax"] > 0.1)]
clsStars = pd.read_csv("./data/clsStars2.csv")

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

occurrenceRates = np.zeros((4,4,6,3))
rawOccurrence = np.zeros((4,4,6,3))
clsOccurrenceRates = np.zeros((4,4,3))
clsRawOccurrence = np.zeros((4,4,3))

colours = ["tab:purple", "tab:pink", "tab:red", "tab:cyan", "tab:green", "tab:orange"]
massBinLabels = ["FGK", "K", "G", "F"]
metBinLabels = ["All [Fe/H]", "[Fe/H] < -0.1", "-0.1 <= [Fe/H] < 0.1", "[Fe/H] >= 0.1"]
#metBinLabels = ["All [Fe/H]", "[Fe/H] < 0.1", "[Fe/H] >= 0.1"]
companionLabels = ["CJ", "WJ", "HJ", "CS", "HS", "SE"]

for i in range(4):
    if i == 0:
        massCut = planets.loc[(planets["st_mass"] > 0.55) & (planets["st_mass"] < 1.6)]
        clsMassCut = clsPlanets.loc[(clsPlanets["st_mass"] > 0.55) & (clsPlanets["st_mass"] < 1.6)]
        clsStarMassCut = clsStars.loc[(clsStars["mass"] > 0.55) & (clsStars["mass"] < 1.6)]
    elif i == 1:
        #massCut = planets.loc[planets["st_mass"] < 0.55]
        massCut = planets.loc[(planets["st_mass"] > 0.55) & (planets["st_mass"] < 0.8)]
        clsMassCut = clsPlanets.loc[(clsPlanets["st_mass"] > 0.55) & (clsPlanets["st_mass"] < 0.8)]
        clsStarMassCut = clsStars.loc[(clsStars["mass"] > 0.55) & (clsStars["mass"] < 0.8)]
    #elif i == 2: 
        #massCut = planets.loc[(planets["st_mass"] >= 0.55) & (planets["st_mass"] < 0.8)]
        #massCut = planets.loc[planets["st_mass"] >= 1]
    elif i == 2:
        massCut = planets.loc[(planets["st_mass"] >= 0.8) & (planets["st_mass"] < 1.05)]
        clsMassCut = clsPlanets.loc[(clsPlanets["st_mass"] >= 0.9) & (clsPlanets["st_mass"] < 1.05)]
        clsStarMassCut = clsStars.loc[(clsStars["mass"] >= 0.9) & (clsStars["mass"] < 1.1)]
    elif i == 3:
        massCut = planets.loc[(planets["st_mass"] >= 1.05) & (planets["st_mass"] < 1.6)]
        clsMassCut = clsPlanets.loc[(clsPlanets["st_mass"] >= 1.05) & (clsPlanets["st_mass"] < 1.6)]
        clsStarMassCut = clsStars.loc[(clsStars["mass"] >= 1.05) & (clsStars["mass"] < 1.6)]
    
    for j in range(4):
        if j == 0:
            metCut = massCut
            clsMetCut = clsMassCut
            clsStarMetCut = clsStarMassCut
        elif j == 1:
            metCut = massCut.loc[massCut["st_met"] < -0.1]
            clsMetCut = clsMassCut.loc[clsMassCut["st_met"] < -0.1]
            clsStarMetCut = clsStarMassCut.loc[clsStarMassCut["[Fe/H]"] < -0.1]
            #metCut = massCut.loc[massCut["st_met"] < 0]
        elif j == 2:
            metCut = massCut.loc[(massCut["st_met"] >= -0.1) & (massCut["st_met"] < 0.1)]
            clsMetCut = clsMassCut.loc[(clsMassCut["st_met"] >= -0.1) & (clsMassCut["st_met"] < 0.1)]
            clsStarMetCut = clsStarMassCut.loc[(clsStarMassCut["[Fe/H]"] >= -0.1) & (clsStarMassCut["[Fe/H]"] < 0.1)]
        #    #metCut = massCut.loc[massCut["st_met"] >= 0]
        elif j == 3:
            metCut = massCut.loc[massCut["st_met"] >= 0.1]
            clsMetCut = clsMassCut.loc[clsMassCut["st_met"] >= 0.1]
            clsStarMetCut = clsStarMassCut.loc[clsStarMassCut["[Fe/H]"] >= 0.1]
        print("CLS", massBinLabels[i], metBinLabels[j])
        clsOccurrenceRates[i,j], clsRawOccurrence[i,j] = occurrenceRate(clsMetCut, clsStarMetCut, minMassIdx, maxMassIdx, minSmaxIdx, maxSmaxIdx, clsFlag=True)
        

        for k in range(6):
            if k == 0:
                hostPopulation = metCut.loc[metCut["pl_type"] == "CJ"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 13 == 0]
            if k == 1:
                hostPopulation = metCut.loc[metCut["pl_type"] == "WJ"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 11 == 0]
            if k == 2:
                hostPopulation = metCut.loc[metCut["pl_type"] == "HJ"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 7 == 0]
            if k == 3:
                hostPopulation = metCut.loc[metCut["pl_type"] == "CS"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 5 == 0]
            if k == 4:
                hostPopulation = metCut.loc[metCut["pl_type"] == "HS"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 3 == 0]
            if k == 5:
                hostPopulation = metCut.loc[metCut["pl_type"] == "SE"]
                targetPopulation = metCut.loc[metCut["companion_type"] % 2 == 0]
            print(companionLabels[k], massBinLabels[i], metBinLabels[j])
            occurrenceRates[i,j,k], rawOccurrence[i,j,k] = occurrenceRate(targetPopulation, hostPopulation, minMassIdx, maxMassIdx, minSmaxIdx, maxSmaxIdx)

fig, ax = plt.subplots(6,4, figsize = (24,24))
for i in range(4):
    #for ii in range(1):
    #    iii = 2*i+ii+1
    for k in range(6):
                #print(iii,j,k,occurrenceRates[iii,j,k,0])
        ax[k][i].bar(metBinLabels, occurrenceRates[i,:,k,0], yerr = np.reshape(np.array([occurrenceRates[i,:,k,1:]]).T, shape = (2,4)), capsize = 5, color = colours[k], alpha = 0.5)
        ax[k][i].bar(metBinLabels, clsOccurrenceRates[i,:,0], yerr = np.reshape(np.array([clsOccurrenceRates[i,:,1:]]).T, shape = (2,4)), capsize = 5, ecolor = "tab:blue", edgecolor = "tab:blue", fill = False)
        ax[k][i].set_ylim(0,0.8)
    ax[0][i].set_title(massBinLabels[i])

plt.show()

fig1, ax1 = plt.subplots(1,1)

im = ax1.imshow(clsOccurrenceRates[:,:,0], vmin = 0 , vmax = 0.6)
ax1.set_yticks(range(len(massBinLabels)), labels = massBinLabels, rotation = 45)
ax1.set_xticks(range(len(metBinLabels)), labels = metBinLabels)
ax1.set_title("CLS Gas Giant Occurrence Rates")

for i in range(len(massBinLabels)):
    for j in range(len(metBinLabels)):
        ax1.text(j, i-0.2, "{0:.2f}".format(clsOccurrenceRates[i,j,0]), color = "white", ha = "center")
        ax1.text(j, i+0.2, f"{int(clsRawOccurrence[i,j,0])} / {int(clsRawOccurrence[i,j,1])}", color = "white", ha = "center")

fig1.savefig("./plots/clsGiantOccurenceRates.png")
plt.show()

for k in range(6):
    fig2, ax2 = plt.subplots(1,1, figsize = (6,6))

    im = ax2.imshow(occurrenceRates[:,:,k,0], vmin = 0 , vmax = 0.5)
    ax2.set_yticks(range(len(massBinLabels)), labels = massBinLabels, rotation = 45)
    ax2.set_xticks(range(len(metBinLabels)), labels = metBinLabels)
    ax2.set_title(companionLabels[k]+" Occurrence Rates")

    for i in range(len(massBinLabels)):
        for j in range(len(metBinLabels)):
            ax2.text(j, i-0.2, "{0:.2f}".format(occurrenceRates[i,j,k,0]), color = "white", ha = "center")
            ax2.text(j, i+0.2, f"{int(rawOccurrence[i,j,k,0])} / {int(rawOccurrence[i,j,k,1])} - {int(10*rawOccurrence[i,j,k,2])/10}", color = "white", ha = "center")

    plt.show()