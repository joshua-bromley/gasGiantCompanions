import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats as stats
from scipy import special as sp

def occurrence(x, a, b):
    return (1/sp.beta(a,b))*x**(a-1) * (1-x)**(b-1)

clsStars = pd.read_csv("./data/clsStars.csv")
clsStars = clsStars.loc[clsStars["mass"] >= 0.6]
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmasse"] > 0.5*317) & (clsPlanets["pl_orbsmax"] > 1) & (clsPlanets["st_mass"] >= 0.6)]



planets = pd.read_csv("./data/gasGiantData.csv")
planets = planets.loc[(planets["st_mass"] > 0.6)]

superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbsmax"] > 1) & (planets["pl_bmasse"] > 0.5*317)]

outerCompanions = planets.loc[planets["companion_type"] > 1] 

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
ggCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]


metallicityCutoffs = np.arange(-0.2,0.31, 0.01)
spotlightPoints = [20,22,24]


SEmetalRichEnhancement = np.zeros_like(metallicityCutoffs)
SEmetalPoorEnhancement = np.zeros_like(metallicityCutoffs)


CompMetalRichEnhancement = np.zeros_like(metallicityCutoffs)
CompMetalPoorEnhancement = np.zeros_like(metallicityCutoffs)

GGMetalRichEnhancement = np.zeros_like(metallicityCutoffs)
GGMetalPoorEnhancement = np.zeros_like(metallicityCutoffs)

clsMetMedian = np.median(clsStars["[Fe/H]"].dropna())

print(f"CLS Median: {clsMetMedian}")
planetsMetMedian = np.median(superEarths.drop_duplicates(subset="hostname")["st_met"].dropna())
print(f"B&B Median: {planetsMetMedian}")


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
    occurRateSEMetalRich = stats.beta.median(nMetalRichSECompanions + 1, nMetalRichSE - nMetalRichSECompanions + 1)
    occurRateSEMetalPoor = stats.beta.median(nMetalPoorSECompanions + 1, nMetalPoorSE - nMetalPoorSECompanions + 1)
    occurRateCompMetalRich = stats.beta.median(nMetalRichCompanions + 1, nMetalRichInner - nMetalRichCompanions + 1)
    occurRateCompMetalPoor = stats.beta.median(nMetalPoorCompanions + 1, nMetalPoorInner - nMetalPoorCompanions +1)
    occurRateGGMetalRich = stats.beta.median(nMetalRichGGCompanions + 1, nMetalRichGG - nMetalRichGGCompanions + 1)
    occurRateGGMetalPoor = stats.beta.median(nMetalPoorGGCompanions + 1, nMetalPoorGG - nMetalPoorGGCompanions + 1)
    
    if occurRateSEMetalRich > occurRateCLSMetalRich:
        sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1] - occurRateCLSMetalRich
        sigmaSEMetalRich = occurRateSEMetalRich - stats.beta.interval(0.68, nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)[0]
    else:
        sigmaCLSMetalRich = occurRateCLSMetalRich - stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[0]
        sigmaSEMetalRich = stats.beta.interval(0.68, nMetalRichSECompanions, nMetalRichSE - nMetalRichSECompanions + 1)[1] - occurRateSEMetalRich
    
    if occurRateSEMetalPoor > occurRateCLSMetalPoor:
        sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)[1] - occurRateCLSMetalPoor
        sigmaSEMetalPoor= occurRateSEMetalPoor - stats.beta.interval(0.68, nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)[0]
    else:
        sigmaCLSMetalPoor = occurRateCLSMetalPoor -  stats.beta.interval(0.68,nclsMetalPoorStars + 1, nclsMetalPoorPlanets - nclsMetalPoorPlanets + 1)[0]
        sigmaSEMetalPoor= stats.beta.interval(0.68, nMetalPoorSECompanions, nMetalPoorSE - nMetalPoorSECompanions + 1)[1] - occurRateSEMetalPoor

    if occurRateCompMetalRich > occurRateCLSMetalRich:
        sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1] - occurRateCLSMetalRich
        sigmaCompMetalRich = occurRateCompMetalRich - stats.beta.interval(0.68, nMetalRichCompanions, nMetalRichInner - nMetalRichCompanions + 1)[0]
    else:
        sigmaCLSMetalRich = occurRateCLSMetalRich - stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[0]
        sigmaCompMetalRich = stats.beta.interval(0.68, nMetalRichCompanions, nMetalRichInner - nMetalRichCompanions + 1)[1] - occurRateCompMetalRich
    
    if occurRateCompMetalPoor > occurRateCLSMetalPoor:
        sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)[1] - occurRateCLSMetalPoor
        sigmaCompMetalPoor= occurRateCompMetalPoor - stats.beta.interval(0.68, nMetalPoorCompanions, nMetalPoorInner - nMetalPoorCompanions + 1)[0]
    else:
        sigmaCLSMetalPoor = occurRateCLSMetalPoor -  stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)[0]
        sigmaCompMetalPoor= stats.beta.interval(0.68, nMetalPoorCompanions, nMetalPoorInner - nMetalPoorCompanions + 1)[1] - occurRateCompMetalPoor

    if occurRateGGMetalRich > occurRateCLSMetalRich:
        sigmaCLSMetalRich = stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[1] - occurRateCLSMetalRich
        sigmaGGMetalRich = occurRateGGMetalRich - stats.beta.interval(0.68, nMetalRichGGCompanions, nMetalRichGG - nMetalRichGGCompanions + 1)[0]
    else:
        sigmaCLSMetalRich = occurRateCLSMetalRich - stats.beta.interval(0.68,nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1)[0]
        sigmaGGMetalRich = stats.beta.interval(0.68, nMetalRichGGCompanions, nMetalRichGG - nMetalRichGGCompanions + 1)[1] - occurRateGGMetalRich
    
    if occurRateGGMetalPoor > occurRateCLSMetalPoor:
        sigmaCLSMetalPoor = stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)[1] - occurRateCLSMetalPoor
        sigmaGGMetalPoor= occurRateGGMetalPoor - stats.beta.interval(0.68, nMetalPoorGGCompanions, nMetalPoorGG - nMetalPoorGGCompanions + 1)[0]
    else:
        sigmaCLSMetalPoor = occurRateCLSMetalPoor -  stats.beta.interval(0.68,nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1)[0]
        sigmaGGMetalPoor= stats.beta.interval(0.68, nMetalPoorGGCompanions, nMetalPoorGG - nMetalPoorGGCompanions + 1)[1] - occurRateGGMetalPoor
    


        
    SEmetalRichEnhancement[i] = (occurRateSEMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaSEMetalRich*sigmaSEMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    SEmetalPoorEnhancement[i] = (occurRateSEMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaSEMetalPoor*sigmaSEMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)

    CompMetalRichEnhancement[i] = (occurRateCompMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaCompMetalRich*sigmaCompMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    CompMetalPoorEnhancement[i] = (occurRateCompMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaCompMetalPoor*sigmaCompMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)

    GGMetalRichEnhancement[i] = (occurRateGGMetalRich - occurRateCLSMetalRich)/np.sqrt(sigmaGGMetalRich*sigmaGGMetalRich + sigmaCLSMetalRich*sigmaCLSMetalRich)
    GGMetalPoorEnhancement[i] = (occurRateGGMetalPoor - occurRateCLSMetalPoor)/np.sqrt(sigmaGGMetalPoor*sigmaGGMetalPoor + sigmaCLSMetalPoor*sigmaCLSMetalPoor)

   

    if i in spotlightPoints:
        if i == 22:
            cutoff = 0.023
        elif i == 24:
            cutoff = 0.041
        x = np.linspace(0,0.5,1000)
        #fig1, ax1 = plt.subplots(1,1)
        #ax1.plot(x, occurrence(x, nclsMetalRichPlanets + 1, nclsMetalRichStars - nclsMetalRichPlanets + 1), ls = "-", color = "tab:blue", label = f"P(GG|[Fe/H] > {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nclsMetalPoorPlanets + 1, nclsMetalPoorStars - nclsMetalPoorPlanets + 1), ls = "-", color = "tab:orange", label = f"P(GG|[Fe/H] <= {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nMetalRichSECompanions + 1, nMetalRichSE - nMetalRichSECompanions + 1), ls = "--", color = "tab:blue", label = f"P(GG|SE,[Fe/H] > {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nMetalPoorSECompanions + 1, nMetalPoorSE - nMetalPoorSECompanions + 1), ls = "--", color = "tab:orange", label = f"P(GG|SE,[Fe/H] <= {np.round(cutoff, decimals=1)})")
        #ax1.legend(frameon = False)
        plt.tight_layout()
        #fig1.savefig(f"./plots/occurRateMetCutoff{np.round(cutoff, decimals=2)}.png")
        print(f"Cutoff: {np.round(cutoff, decimals=5)}")
        print(f"Metal Rich Field: {nclsMetalRichPlanets}/{nclsMetalRichStars}, {occurRateCLSMetalRich} +/- {sigmaCLSMetalRich}")
        print(f"Metal Rich GG: {nMetalRichGGCompanions}/{nMetalRichGG}, {occurRateGGMetalRich} +/- {sigmaGGMetalRich}")
        print(f"Metal Poor Field: {nclsMetalPoorPlanets}/{nclsMetalPoorStars}, {occurRateCLSMetalPoor} +/- {sigmaCLSMetalPoor}")
        print(f"Metal Poor GG: {nMetalPoorGGCompanions}/{nMetalPoorGG}, {occurRateGGMetalPoor} +/- {sigmaGGMetalPoor}")
        print(f"Metal Rich Enhancement: {GGMetalRichEnhancement[i]}")
        print(f"Metal Poor Enhancement: {GGMetalPoorEnhancement[i]}")



 


fig, ax = plt.subplots(1,1)
ax.plot(metallicityCutoffs, SEmetalRichEnhancement, ls = "", marker = "o", label = "[Fe/H] > x, SE Inner")
ax.plot(metallicityCutoffs, SEmetalPoorEnhancement, ls = "", marker = "o", label = "[Fe/H] <= x, SE Inner")
ax.plot(metallicityCutoffs, GGMetalRichEnhancement, ls = "", marker = "o", label = "[Fe/H] > x, Any Inner")
ax.plot(metallicityCutoffs, GGMetalPoorEnhancement, ls = "", marker = "o", label = "[Fe/H] <= x, Any Inner")
ax.vlines([clsMetMedian, planetsMetMedian], 0,5)
ax.legend(frameon = False)
ax.set_ylabel("Enhancement ($\sigma$)")
ax.set_xlabel('Metallicity Cutoff ([Fe/H])')
plt.show()

#fig.savefig("./plots/metallicityEnhancement.png")

#print(metallicityCutoffs)

    




