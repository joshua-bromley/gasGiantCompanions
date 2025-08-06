import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import stats as stats

planets = pd.read_csv("./data/gasGiantData.csv")

hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0) & (pd.isna(planets["pl_orbsmax"]) == False)]

hjCompanionMass = np.log10(np.array(hjCompanions["pl_bmassj"].values))
hjCompanionSMax = np.log10(np.array(hjCompanions["pl_orbsmax"].values))

hjCompanionLogData = np.vstack((hjCompanionMass, hjCompanionSMax))

kernel = stats.gaussian_kde(hjCompanionLogData)

logMassBins = np.linspace(np.log10(0.5), np.log10(20))
logSmaxBins = np.linspace(-1,np.log10(20))
x,y = np.meshgrid(logMassBins, logSmaxBins)
positions = np.vstack([x.ravel(), y.ravel()])
z = np.reshape(kernel(positions), x.shape)
z = z/np.max(z)


wcjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)] 
wcjCompanions = wcjCompanions.loc[(pd.isna(wcjCompanions["pl_orbsmax"]) == False)]
wcjLogMass = np.log10(np.array(wcjCompanions["pl_bmassj"]))
wcjLogSmax = np.log10(np.array(wcjCompanions["pl_orbsmax"]))



correctWCJSample = pd.DataFrame(columns=planets.columns)

def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return array[idx]

for i in range(len(wcjCompanions)):
    planetMass = np.log10(wcjCompanions.iloc[i]["pl_bmassj"])
    planetSmax = np.log10(wcjCompanions.iloc[i]["pl_orbsmax"])
    massIdx = np.abs(logMassBins - planetMass).argmin()
    smaxIdx = np.abs(logSmaxBins - planetSmax).argmin()
    probability = z[massIdx, smaxIdx]
    randNum = np.random.rand(1)
    if probability > randNum:
        correctWCJSample.loc[len(correctWCJSample)] = wcjCompanions.iloc[i]

print(len(correctWCJSample))
correctWCJLogMass = np.log10(np.array(correctWCJSample["pl_bmassj"]).astype(float))
correctWCJLogSmax = np.log10(np.array(correctWCJSample["pl_orbsmax"]).astype(float))

coldJupiters = planets.loc[(planets["pl_type"] == "CJ") & (pd.isna(planets["pl_orbsmax"]) == False)]
cjLogMass = np.log10(np.array(coldJupiters["pl_bmassj"]))
cjLogSmax = np.log10(np.array(coldJupiters["pl_orbsmax"]))

correctCJSample = pd.DataFrame(columns=planets.columns)

for i in range(len(coldJupiters)):
    planetMass = np.log10(coldJupiters.iloc[i]["pl_bmassj"])
    planetSmax = np.log10(coldJupiters.iloc[i]["pl_orbsmax"])
    massIdx = np.abs(logMassBins - planetMass).argmin()
    smaxIdx = np.abs(logSmaxBins - planetSmax).argmin()
    probability = z[massIdx, smaxIdx]
    randNum = np.random.rand(1)
    if probability > randNum:
        correctCJSample.loc[len(correctCJSample)] = coldJupiters.iloc[i]

print(len(correctCJSample))
correctCJLogMass = np.log10(np.array(correctCJSample["pl_bmassj"]).astype(float))
correctCJLogSmax = np.log10(np.array(correctCJSample["pl_orbsmax"]).astype(float))

fig, ax = plt.subplots(1,1)
ax.imshow(z, origin = "lower", extent = [np.log10(0.5), np.log10(20), -1, np.log10(20)], vmax = 1, cmap = "cividis")
ax.plot(hjCompanionMass, hjCompanionSMax, ls ="", marker = "o", color = "tab:red")
ax.plot(correctWCJLogMass,correctWCJLogSmax , ls = "", marker = "o", color = "tab:orange")
ax.plot(wcjLogMass, wcjLogSmax, ls = "", marker = "+", color = "tab:green")
ax.plot(correctCJLogMass, correctCJLogSmax, ls = "", marker = ".", color = "tab:pink")
plt.show()

fig, ax = plt.subplots(1,1)
bins = np.arange(0.0,1.001,0.1)
ax.hist(hjCompanions["pl_orbeccen"], bins = bins, histtype="step",density = False, label = "HJ Companions", color = "tab:red")
ax.hist(correctWCJSample["pl_orbeccen"], bins = bins, histtype = "step",density = False, label = "Matched WCJ Companions", color = "tab:orange")
#ax.hist(wcjCompanions["pl_orbeccen"], bins = bins, density=True, histtype="step", label ="WCJ Companions", color = "tab:green")
ax.hist(correctCJSample["pl_orbeccen"], bins = bins, histtype = "step",density = False, label = "Matched Cold Jupiters", color = "tab:pink")

ax.legend(frameon = False)
plt.show()