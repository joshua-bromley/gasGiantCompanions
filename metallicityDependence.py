import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats as stats

clsStars = pd.read_csv("./data/clsStars.csv")
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmasse"] > 0.5*317) & (clsPlanets["pl_orbsmax"] > 1)]

planets = pd.read_csv("./data/gasGiantData.csv")

superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbsmax"] > 1)]

outerCompanions = planets.loc[planets["companion_type"] > 1] 

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
ggCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]


metallicityCutoffs = np.linspace(-0.2,0.3)


SEmetalRichEnhancement = np.zeros_like(metallicityCutoffs)
SEmetalPoorEnhancement = np.zeros_like(metallicityCutoffs)


CompMetalRichEnhancement = np.zeros_like(metallicityCutoffs)
CompMetalPoorEnhancement = np.zeros_like(metallicityCutoffs)

GGmetalRichEnhancement = np.zeros_like(metallicityCutoffs)
GGmetalPoorEnhancement = np.zeros_like(metallicityCutoffs)

clsMetMedian = np.median(clsStars["[Fe/H]"].dropna())

print(np.min(clsStars["[Fe/H]"]), np.max(clsStars["[Fe/H]"]))


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

    nMetalRichCompanions = len(np.unique(outerCompanions.loc[outerCompanions["st_met"] > cutoff]["hostname"]))
    nMetalPoorCompanions = len(np.unique(outerCompanions.loc[outerCompanions["st_met"] <= cutoff]["hostname"]))
    nMetalRichInner = len(np.unique(planets.loc[planets["st_met"] > cutoff]["hostname"]))
    nMetalPoorInner = len(np.unique(planets.loc[planets["st_met"] <= cutoff]["hostname"]))

    nMetalRichGG = len(np.unique(gasGiants.loc[gasGiants["st_met"] > cutoff]["hostname"]))
    nMetalPoorGG = len(np.unique(gasGiants.loc[gasGiants["st_met"] <= cutoff]["hostname"]))
    nMetalRichGGCompanions = len(np.unique(ggCompanions.loc[ggCompanions["st_met"] > cutoff]["hostname"]))
    nMetalPoorGGCompanions = len(np.unique(ggCompanions.loc[ggCompanions["st_met"] <= cutoff]["hostname"]))

    occurRateCLSMetalRich = stats.beta.median(nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)
    occurRateCLSMetalPoor = stats.beta.median(nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)
    occurRateSEMetalRich = stats.beta.median(nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)
    occurRateSEMetalPoor = stats.beta.median(nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)
    occurRateCompMetalRich = stats.beta.median(nMetalRichCompanions + 1, nMetalRichInner - nMetalRichCompanions + 1)
    occurRateCompMetalPoor = stats.beta.median(nMetalPoorCompanions + 1, nMetalPoorInner - nMetalPoorCompanions +1)
    occurRateGGMetalRich = stats.beta.median(nMetalRichGGCompanions, nMetalRichGG - nMetalRichGGCompanions + 1)
    occurRateGGMetalPoor = stats.beta.median(nMetalPoorGGCompanions, nMetalPoorGG - nMetalPoorGGCompanions + 1)
    

    if occurRateSEMetalRich > occurRateCLSMetalRich:
        sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1] - occurRateCLSMetalRich
        sigmaSEMetalRich = occurRateSEMetalRich - stats.beta.interval(0.68, nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)[0]
    else:
        sigmaCLSMetalRich = occurRateCLSMetalRich - stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[0]
        sigmaSEMetalRich = stats.beta.interval(0.68, nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)[1] - occurRateSEMetalRich
    
    if occurRateSEMetalPoor > occurRateCLSMetalPoor:
        sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[1] - occurRateCLSMetalPoor
        sigmaSEMetalPoor= occurRateSEMetalPoor - stats.beta.interval(0.68, nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)[0]
    else:
        sigmaCLSMetalPoor = occurRateCLSMetalPoor -  stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[0]
        sigmaSEMetalPoor= stats.beta.interval(0.68, nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)[1] - occurRateSEMetalPoor

    if occurRateCompMetalRich > occurRateCLSMetalRich:
        sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1] - occurRateCLSMetalRich
        sigmaCompMetalRich = occurRateCompMetalRich - stats.beta.interval(0.68, nMetalRichCompanions, nMetalRichInner - nMetalRichCompanions + 1)[0]
    else:
        sigmaCLSMetalRich = occurRateCLSMetalRich - stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[0]
        sigmaCompMetalRich = stats.beta.interval(0.68, nMetalRichCompanions, nMetalRichInner - nMetalRichCompanions + 1)[1] - occurRateCompMetalRich
    
    if occurRateCompMetalPoor > occurRateCLSMetalPoor:
        sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[1] - occurRateCLSMetalPoor
        sigmaCompMetalPoor= occurRateCompMetalPoor - stats.beta.interval(0.68, nMetalPoorCompanions, nMetalPoorInner - nMetalPoorCompanions + 1)[0]
    else:
        sigmaCLSMetalPoor = occurRateCLSMetalPoor -  stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[0]
        sigmaCompMetalPoor= stats.beta.interval(0.68, nMetalPoorCompanions, nMetalPoorInner - nMetalPoorCompanions + 1)[1] - occurRateCompMetalPoor


        
    SEmetalRichEnhancement[i] = (occurRateSEMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaSEMetalRich*sigmaSEMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    SEmetalPoorEnhancement[i] = (occurRateSEMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaSEMetalPoor*sigmaSEMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)

    CompMetalRichEnhancement[i] = (occurRateCompMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaCompMetalRich*sigmaCompMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    CompMetalPoorEnhancement[i] = (occurRateCompMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaCompMetalPoor*sigmaCompMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)
 


fig, ax = plt.subplots(1,1)
ax.plot(metallicityCutoffs, SEmetalRichEnhancement, ls = "", marker = "o", label = "[Fe/H] > x, SE Inner")
ax.plot(metallicityCutoffs, SEmetalPoorEnhancement, ls = "", marker = "o", label = "[Fe/H] <= x, SE Inner")
#ax.plot(metallicityCutoffs, CompMetalRichEnhancement, ls = "", marker = "o", label = "[Fe/H] > x, Any Inner")
#ax.plot(metallicityCutoffs, CompMetalPoorEnhancement, ls = "", marker = "o", label = "[Fe/H] <= x, Any Inner")
ax.vlines(clsMetMedian, 0,2.5)
ax.legend(frameon = False)
ax.set_ylabel("Enhancement ($\sigma$)")
ax.set_xlabel('Metallicity Cutoff ([Fe/H])')

fig.savefig("./plots/metallicityEnhancement.png")

    




