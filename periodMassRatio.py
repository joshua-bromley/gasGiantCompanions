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

planets = pd.read_csv("./data/gasGiantDataComplete.csv").drop_duplicates(subset = ["pl_name"])
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]


planetPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
hjPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
wcjPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
sePairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
ssPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])


hostnames = np.unique(planets["hostname"])

for i in range(len(hostnames)):
    hostname = hostnames[i]
    system = planets.loc[planets["hostname"] == hostname].sort_values(by="pl_orbper")
    for j in range(1,len(system)):
        planetA = system.iloc[j]
        planetB = system.iloc[j-1]
        planetAName = planetA["pl_name"]
        planetBName = planetB["pl_name"]
        planetAMass = planetA["pl_bmasse"]
        planetAPer = planetA["pl_orbper"]
        planetAeccen = planetA["pl_orbeccen"]
        planetBMass = planetB["pl_bmasse"]
        planetBPer = planetB["pl_orbper"]
        planetBeccen = planetB["pl_orbeccen"]
        stellarMass = planetA["st_mass"]
        stellarMet = planetA["st_met"]
        massRatio = planetAMass/planetBMass
        perRatio = planetAPer/planetBPer
        planetPairs.loc[len(planetPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
        if planetB["pl_type"] == "HJ" and (planetA["pl_type"] == "WJ" or planetA["pl_type"] == "CJ"):
            hjPairs.loc[len(hjPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
        elif (planetB["pl_type"] == "WJ" or planetB["pl_type"] == "CJ") and (planetA["pl_type"] == "WJ" or planetA["pl_type"] == "CJ"):
            wcjPairs.loc[len(wcjPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
        elif (planetB["pl_type"] == "SE") and (planetA["pl_type"] == "CJ" or planetA["pl_type"] == "WJ"):
            sePairs.loc[len(sePairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
        elif (planetB["pl_type"] == "HS" or planetB["pl_type"] == "CS") and (planetA["pl_type"] == "WJ" or planetA["pl_type"] == 'CJ'):
            ssPairs.loc[len(ssPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]


 

highMassPairs = planetPairs.loc[planetPairs["pla_mass"] > 0.5*317]

periodRatios = planetPairs["per_ratio"].values
massRatios = planetPairs["mass_ratio"].values
outerEccen = planetPairs["pla_eccen"].values

hjPerRatios = hjPairs["per_ratio"].values
hjMassRatios = hjPairs["mass_ratio"].values

wcjPerRatios = wcjPairs["per_ratio"].values
wcjMassRatios = wcjPairs["mass_ratio"].values

perBins = [1, 3, 10, 30, 100, 300, 1000, 3000]
perBinCentres = np.sqrt(np.multiply(perBins[:-1], perBins[1:]))
massBins = [0.1, 0.3, 1, 3, 10, 30]
massBinCentres = np.sqrt(np.multiply(massBins[:-1], massBins[1:]))

perRatioDist = np.zeros((4, len(perBinCentres), 3))
massRatioDist = np.zeros((2,len(massBinCentres), 3))

planetPairsList = [sePairs, ssPairs, hjPairs, wcjPairs]
planetPairLabels= ['SE-CJ Pair', "SS-CJ Pair",  "HJ-CJ Pair", "CJ-CJ Pair"]
colours = ["tab:orange", "tab:green", "tab:red", "tab:purple"]
symbols = [ "s", 'h', '*', "P"]

for i in range(len(planetPairsList)):
    pairsList = planetPairsList[i]
    nPairs = len(pairsList)
    perRatioCounts = np.zeros(len(perBinCentres))
    for j in range(len(pairsList)):
        perRatio = pairsList.iloc[j]["per_ratio"]
        for k in range(len(perBinCentres)):
            if perRatio > perBins[k] and perRatio < perBins[k+1]:
                perRatioCounts[k] += 1
    
    for j in range(len(perRatioCounts)):
        a = perRatioCounts[j] + 1
        b = nPairs - perRatioCounts[j] + 1
        perRatioDist[i][j][0] = stats.beta.median(a,b)
        errorBars = stats.beta.interval(0.68, a, b)
        perRatioDist[i][j][1] = perRatioDist[i][j][0] - errorBars[0]
        perRatioDist[i][j][2] = errorBars[1] - perRatioDist[i][j][0]

    if i > 1:
        pairsList = planetPairsList[i]
        nPairs = len(pairsList)
        massRatioCounts = np.zeros(len(massBinCentres))
        for j in range(len(pairsList)):
            massRatio = pairsList.iloc[j]["mass_ratio"]
            for k in range(len(massBinCentres)):
                if massRatio > massBins[k] and massRatio < massBins[k+1]:
                    massRatioCounts[k] += 1
        
        for j in range(len(massRatioCounts)):
            a = massRatioCounts[j] + 1
            b = nPairs - massRatioCounts[j] + 1
            massRatioDist[i-2][j][0] = stats.beta.median(a,b)
            errorBars = stats.beta.interval(0.68, a, b)
            massRatioDist[i-2][j][1] = massRatioDist[i-2][j][0] - errorBars[0]
            massRatioDist[i-2][j][2] = errorBars[1] - massRatioDist[i-2][j][0]


fig, ax = plt.subplots(1,2, figsize = (12, 5))
#ax[0].hist(hjPairs["mass_ratio"], histtype ="step", color = "tab:red", bins = massBins)
#ax[0].hist(wcjPairs["mass_ratio"], histtype = "step", color = "tab:purple", bins = massBins)

#ax[1].hist(hjPairs["per_ratio"], histtype = "step", color = "tab:red", bins = perBins, label = "HJ Inner")
#ax[1].hist(wcjPairs["per_ratio"], histtype = "step", color = "tab:purple", bins = perBins, label = "CJ Inner")
#ax[1].hist(sePairs["per_ratio"], histtype = "step", color = "tab:orange", bins = perBins, label = "SE Inner")
#ax[1].hist(ssPairs["per_ratio"], histtype = "step", color = "tab:green", bins = perBins, label = "SS Inner")

for i in range(len(perRatioDist)):
    ax[1].errorbar(perBinCentres, perRatioDist[i][:,0], yerr = perRatioDist[i][:,1:].T, ls = "", marker = symbols[i], label = planetPairLabels[i], color = colours[i], alpha = 0.8, capsize = 5)

for i in [2,3]:
    ax[0].errorbar(massBinCentres,massRatioDist[i-2][:,0], yerr = massRatioDist[i-2][:,1:].T, ls = "", marker = symbols[i], label = planetPairLabels[i], color = colours[i], alpha = 0.8, capsize = 5)


#for i in range(len(hjPairs)):
    #ax.text(hjPairs["per_ratio"].values[i], hjPairs["mass_ratio"].values[i], hjPairs.iloc[i]["hostname"])
ax[0].set_xscale("log")
#ax.set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_xlabel("Period Ratio", fontsize = 16)
ax[0].set_ylabel("Occurrence", fontsize = 16)
ax[0].set_xlabel('Mass Ratio', fontsize = 16)
ax[1].legend(frameon = False)

yTicks = np.arange(0.0,0.81,0.2)
ax[0].set_yticks(yTicks)
ax[0].set_xlim(0.1,20)
ax[0].set_ylim(0,0.8)

yTicks = np.arange(0.0,0.61,0.1)
ax[1].set_yticks(yTicks)
ax[1].set_xlim(1,3000)
ax[1].set_ylim(0,0.6)

tickLabelSize = 12
ax[0].tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[0].tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax[1].tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()



fig.savefig("./plots/compltCorrPerMassRatio.png")
fig.savefig("./plots/compltCorrPerMassRatio.pdf")



print(np.mean(hjMassRatios), np.std(hjMassRatios)/np.sqrt(len(hjMassRatios)))
print(np.mean(wcjMassRatios), np.std(wcjMassRatios)/np.sqrt(len(wcjMassRatios)))

print(np.mean(hjPerRatios), np.std(hjPerRatios)/np.sqrt(len(hjPerRatios)))
print(np.mean(wcjPerRatios), np.std(wcjPerRatios)/np.sqrt(len(wcjPerRatios)))
print(np.mean(sePairs["per_ratio"]), np.std(sePairs["per_ratio"])/np.sqrt(len(sePairs["per_ratio"])))
print(np.mean(ssPairs["per_ratio"]), np.std(ssPairs["per_ratio"])/np.sqrt(len(ssPairs["per_ratio"])))

        

plt.show()