import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp
from scipy import stats as stats

clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows=101)
clsStars = pd.read_csv("./data/clsStars.csv")
planets = pd.read_csv("./data/gasGiantData.csv")

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

planetPopulations = [clsGasGiants, seCompanions, ssCompanions, hjCompanions, wcjCompanions]
hostPopulations = [clsStars, superEarths, subSaturns, hotJupiters, gasGiants]
labels = ["P(GG)", "P(GG|SE)", "P(GG|SS)", "P(GG|HJ)", "P(GG|WCJ)"]

occurrenceRates = np.zeros((len(planetPopulations),3))
enhancement = np.zeros(len(planetPopulations))

fig, ax = plt.subplots(1,1, figsize = (6,4))

x = np.linspace(0,0.5,1000)

for i in range(len(planetPopulations)):
    n_det = len(np.unique(planetPopulations[i]["hostname"]))
    if i > 0:
        n_tot = len(np.unique(hostPopulations[i]["hostname"]))
    else:
        n_tot = len(np.unique(hostPopulations[i]["CPS identifier"]))
    
    a = n_det -1
    b = n_tot - n_det -1 
    occurrenceRates[i][0] = stats.beta.median(a,b)
    errorBars = stats.beta.interval(0.68, a, b)
    occurrenceRates[i][1] = occurrenceRates[i][0] - errorBars[0]
    occurrenceRates[i][2] = errorBars[1] - occurrenceRates[i][0]
    ax.plot(x, occurrence(x,a,b), label = labels[i])
    if i > 0:
        enhancement[i] = (occurrenceRates[i][0] - occurrenceRates[0][0])/(np.sqrt(occurrenceRates[i][1]**2 + occurrenceRates[0][2]**2))

    print(labels[i])
    print(f"{n_det}/{n_tot}, {occurrenceRates[i][0]} (-{occurrenceRates[i][1]} +{occurrenceRates[i][2]}), {enhancement[i]}")

ax.legend(frameon = False)
ax.set_xlabel("P(GG)")
ax.set_ylabel("PDF")
fig.savefig("./plots/outerCompanionOccurrence.png")

plt.show()


