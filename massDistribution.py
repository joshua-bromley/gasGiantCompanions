import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import stats as stats
from scipy import optimize as opt
import pickle

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

nPlanets = len(planets)

coldJupiters = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
hotJupiters = planets.loc[planets["pl_type"] == "HJ"]
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
cjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]

gasGiantsList = [ hjCompanions, cjCompanions]
hostSystemsList = [hotJupiters, coldJupiters]
gasGiantLabels= [ "HJ Companions", "CJ Companions"]
symbols = [ '*', "P"]
colours = [ "tab:red", "tab:purple"]
#gasGiantsList = [hjCompanions]
cmpltMapMassBins = np.logspace(np.log10(0.3), np.log10(30), 50)
cmpltMapSmaxBins = np.logspace(np.log10(0.3), np.log10(30), 50)


massBins = np.logspace(np.log10(0.5), np.log10(25), 6)
massBinsIdx = np.zeros(len(massBins), dtype=np.int32)
massBinCentres = np.sqrt(np.multiply(massBins[:-1],massBins[1:]))
logMassBinWidth = np.mean(np.log10(np.divide(massBins[1:],massBins[:-1])))
smaxBins = np.logspace(np.log10(0.5), np.log10(20), 6)
smaxBinsIdx = np.zeros(len(smaxBins), dtype=np.int32)
smaxBinCentres = np.sqrt(np.multiply(smaxBins[:-1],smaxBins[1:]))
logSmaxBinWidth = np.mean(np.log10(np.divide(smaxBins[1:],smaxBins[:-1])))

for i in range(len(massBins)):
    massBinsIdx[i] = int(np.abs(massBins[i] - cmpltMapMassBins).argmin())
    smaxBinsIdx[i] = int(np.abs(smaxBins[i] - cmpltMapSmaxBins).argmin())

maxSmaxIdx = int(np.abs(20 - cmpltMapSmaxBins).argmin())

minMassIdx = int(np.abs(0.5 - cmpltMapMassBins).argmin())
maxMassIdx = int(np.abs(25 - cmpltMapMassBins).argmin())

massDist = np.zeros((len(gasGiantsList), len(massBins) - 1 , 3))
unCorrMassDist = np.zeros((len(gasGiantsList),len(massBins) - 1, 3))

smaxDist = np.zeros((len(gasGiantsList), len(smaxBins) - 1 , 3))
unCorrSmaxDist = np.zeros((len(gasGiantsList),len(smaxBins) - 1, 3))


for k in range(len(gasGiantsList)):
    gasGiants = gasGiantsList[k]
    hosts = hostSystemsList[k]
    nPlanets = len(hosts)
    completenessMaps = np.zeros((len(hosts),49,49))
    for i in range(len(hosts)):
        hostname = hosts.iloc[i]["hostname"]
        compProb = pickle.load(open("./completenessMaps/"+hostname+".p", "rb"))
        probs = compProb["prob"]
        completenessMaps[i] = probs





    massPlanetCounts = np.zeros(len(massBins) - 1)
    massMissedPlanetCount = np.zeros(len(massBins) - 1)

    smaxPlanetCounts = np.zeros(len(smaxBins) - 1)
    smaxMissedPlanetCount = np.zeros(len(smaxBins) - 1)

    for i in range(len(gasGiants)):
        mass = gasGiants.iloc[i]["pl_bmassj"]
        smax = gasGiants.iloc[i]["pl_orbsmax"]
        
        binIdx = 0
        smaxIdx = 0
        for j in range(len(massBins) - 1):
            
            if mass > massBins[j] and mass < massBins[j+1]:
                binIdx = j
        for j in range(len(smaxBins) - 1):
            
            if smax > smaxBins[j] and smax < smaxBins[j+1]:
                smaxIdx = j
        massPlanetCounts[binIdx] += 1
        smaxPlanetCounts[smaxIdx] += 1

    for i in range(len(hosts)):
        prob = completenessMaps[i]
        for j in range(len(massBins) - 1):
            total = (massBinsIdx[j+1] - massBinsIdx[j])*maxSmaxIdx
            sumProb = np.sum(prob[0:maxSmaxIdx,massBinsIdx[j]:massBinsIdx[j+1]])
            massMissedPlanetCount[j] += (total - sumProb)/total
        for j in range(len(smaxBins) -1):
            total = (smaxBinsIdx[j+1] - smaxBinsIdx[j])*(maxMassIdx - minMassIdx)
            sumProb = np.sum(prob[smaxBinsIdx[j]:smaxBinsIdx[j+1],minMassIdx:maxMassIdx])
            smaxMissedPlanetCount[j] += (total - sumProb)/total




    for i in range(len(massPlanetCounts)):
        a = massPlanetCounts[i] + 1
        b = nPlanets - massPlanetCounts[i] + 1
        bb = nPlanets - massMissedPlanetCount[i] - massPlanetCounts[i] + 1
        massDist[k][i][0] = stats.beta.median(a,bb)
        errorBars = stats.beta.interval(0.68, a, bb)
        massDist[k][i][1] = massDist[k][i][0] - errorBars[0]
        massDist[k][i][2] = errorBars[1] - massDist[k][i][0]

        unCorrMassDist[k][i][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.68, a, b)
        unCorrMassDist[k][i][1] = unCorrMassDist[k][i][0] - errorBars[0]
        unCorrMassDist[k][i][2] = errorBars[1] - unCorrMassDist[k][i][0]
    
    for i in range(len(smaxPlanetCounts)):
        a = smaxPlanetCounts[i] + 1
        b = nPlanets - smaxPlanetCounts[i] + 1
        bb = nPlanets - smaxMissedPlanetCount[i] - smaxPlanetCounts[i] + 1
        smaxDist[k][i][0] = stats.beta.median(a,bb)
        errorBars = stats.beta.interval(0.68, a, bb)
        smaxDist[k][i][1] = smaxDist[k][i][0] - errorBars[0]
        smaxDist[k][i][2] = errorBars[1] - smaxDist[k][i][0]

        unCorrSmaxDist[k][i][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.68, a, b)
        unCorrSmaxDist[k][i][1] = unCorrSmaxDist[k][i][0] - errorBars[0]
        unCorrSmaxDist[k][i][2] = errorBars[1] - unCorrSmaxDist[k][i][0]

def linearModel(x, m, b):
    return m*x + b
powerLawParams = np.zeros((len(gasGiantsList),2))
logMassBinCentres = np.log10(massBinCentres)

for i in range(len(gasGiantsList)):
    param, cov = opt.curve_fit(linearModel, logMassBinCentres, massDist[i][:,0]/logMassBinWidth, sigma = massDist[i][:,1]/logMassBinWidth)
    powerLawParams[i] = param

x = np.logspace(np.log10(0.6), np.log10(20), 1000)

fig, ax = plt.subplots(1,1, figsize = (6,5))

for i in range(len(gasGiantsList)):
    ax.errorbar(massBinCentres, massDist[i][:,0]/logMassBinWidth, yerr = massDist[i][:,1:].T/logMassBinWidth, ls = "", marker = symbols[i], label = gasGiantLabels[i], alpha = 0.8, capsize = 5, color = colours[i])
    #ax.errorbar(massBinCentres, unCorrMassDist[i][:,0], yerr = unCorrMassDist[i][:,1:].T, ls = "", marker = "o")
    ax.plot(x, linearModel(np.log10(x), *powerLawParams[i]), color = colours[i], alpha = 0.5)

#yTicks = [0,0.1,0.2,0.3,0.4]

#ax.set_yticks(yTicks)
#ax.set_xlim(0.6,20)
#ax.set_ylim(0.0,0.4)

ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_xlabel("Mass (MJ)", fontsize = 16)
ax.set_ylabel("Occurrence", fontsize = 16)
ax.legend(frameon = False)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()
#fig.savefig("./plots/compltCorrMassDist.png")
#fig.savefig("./plots/compltCorrMassDist.pdf")
plt.show()

fig1, ax1 = plt.subplots(1,1, figsize = (6,5))
for i in range(len(gasGiantsList)):
    ax1.errorbar(smaxBinCentres, smaxDist[i][:,0]/logSmaxBinWidth, yerr = smaxDist[i][:,1:].T/logSmaxBinWidth, ls = "", marker = symbols[i], label = gasGiantLabels[i], alpha = 0.8, capsize = 5, color = colours[i])
    #ax1.errorbar(smaxBinCentres, unCorrSmaxDist[i][:,0], yerr = unCorrSmaxDist[i][:,1:].T, ls = "", marker = "o")


#yTicks = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6, 1.8]

#ax1.set_yticks(yTicks)
#ax1.set_xlim(0.6,20)
#ax1.set_ylim(0.0,1.8)

ax1.set_xscale("log")
#ax1.set_yscale("log")
ax1.set_xlabel("Semi Major Axis (AU)", fontsize = 16)
ax1.set_ylabel("Occurrence", fontsize = 16)
ax1.legend(frameon = False)

tickLabelSize = 12
ax1.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax1.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax1.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax1.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()
#fig1.savefig("./plots/compltCorrSmaxDist.png")
#fig1.savefig("./plots/compltCorrSmaxDist.pdf")
plt.show()
    






    


