import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats as stats

clsStars = pd.read_csv("./data/clsStars.csv")
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)

planets = pd.read_csv("./data/gasGiantData.csv")

superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[planets["companion_type"] % 2 == 0]

metallicityCutoffs = np.linspace(-0.2,0.2)


metalRichEnhancement = np.zeros_like(metallicityCutoffs)
metalPoorEnhancement = np.zeros_like(metallicityCutoffs)

for i in range(len(metallicityCutoffs)):
    cutoff = metallicityCutoffs[i]
    nclsMetalRichStars = len(np.unique(clsStars.loc[clsStars["[Fe/H]"] > cutoff]["CPS identifier"]))
    nclsMetalPoorStars = len(np.unique(clsStars.loc[clsStars["[Fe/H]"] <= cutoff]["CPS identifier"]))

    nclsMetalRichPlanets = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] > cutoff]["hostname"]))
    nclsMetalPoorPlanets = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] <= cutoff]["hostname"]))

    nMetalRichSE = len(np.unique(superEarths.loc[superEarths["st_met"] > cutoff]["hostname"]))
    nMetalPoorSE = len(np.unique(superEarths.loc[superEarths["st_met"] <= cutoff]["hostname"]))

    nMetalRichSECompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] > cutoff]["hostname"]))
    nMetalPoorSECompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] <= cutoff]["hostname"]))

    occurRateCLSMetalRich = stats.beta.median(nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)
    sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1]

    occurRateCLSMetalPoor = stats.beta.median(nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)
    sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[1]

    occurRateSEMetalRich = stats.beta.median(nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)
    sigmaSEMetalRich = stats.beta.interval(0.68, nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)[0]
    
    occurRateSEMetalPoor = stats.beta.median(nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)
    sigmaSEMetalPoor= stats.beta.interval(0.68, nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)[0]
    metalRichEnhancement[i] = (occurRateSEMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaSEMetalRich*sigmaSEMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    metalPoorEnhancement[i] = (occurRateSEMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaSEMetalPoor*sigmaSEMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)


fig, ax = plt.subplots(1,1)
ax.plot(metallicityCutoffs, metalRichEnhancement, ls = "", marker = "o", label = "[Fe/H] > x")
ax.plot(metallicityCutoffs, metalPoorEnhancement, ls = "", marker = "o", label = "[Fe/H] <= x")
ax.legend(frameon = False)
ax.set_ylabel("Enhancement ($\sigma$)")
ax.set_xlabel('Metallicity Cutoff ([Fe/H])')

fig.savefig("./plots/metallicityEnhancement.png")

    




