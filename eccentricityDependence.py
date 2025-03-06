import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats as stats
from scipy import special as sp

def occurrence(x, a, b):
    return (1/sp.beta(a,b))*x**(a-1) * (1-x)**(b-1)

clsStars = pd.read_csv("./data/clsStars.csv")
clsStars = clsStars.loc[(clsStars["mass"] >= 0.6) & (clsStars["[Fe/H]"] > 0)]
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmasse"] > 0.5*317) & (clsPlanets["pl_orbsmax"] > 1) & (clsPlanets["st_mass"] >= 0.6)  & (clsPlanets["st_met"] > 0)]



planets = pd.read_csv("./data/gasGiantData.csv")
planets = planets.loc[(planets["st_mass"] > 0.6) & (planets["st_met"] > 0)]

superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbsmax"] > 1) & (planets["pl_bmasse"] > 0.5*317)]

outerCompanions = planets.loc[planets["companion_type"] > 1] 

gasGiants = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
ggCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]


eccentricityCutoffs = np.arange(0,1.001,0.1)
spotlightPoints = []


SEDynHotEnhancement = np.zeros_like(eccentricityCutoffs)
SEDynColdEnhancement = np.zeros_like(eccentricityCutoffs)


CompDynHotEnhancement = np.zeros_like(eccentricityCutoffs)
CompDynColdEnhancement = np.zeros_like(eccentricityCutoffs)

GGDynHotEnhancement = np.zeros_like(eccentricityCutoffs)
GGDynColdEnhancement = np.zeros_like(eccentricityCutoffs)

clsMetMedian = np.median(clsStars["[Fe/H]"].dropna())


planetsMetMedian = np.median(superEarths.drop_duplicates(subset="hostname")["st_met"].dropna())



for i in range(len(eccentricityCutoffs)):
    cutoff = eccentricityCutoffs[i]
    nclsDynHotStars = len(np.unique(clsStars.loc[clsStars["[Fe/H]"] > cutoff]["CPS identifier"]))
    nclsDynColdStars = len(np.unique(clsStars.loc[clsStars["[Fe/H]"] <= cutoff]["CPS identifier"]))

    nclsDynHotPlanets = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] > cutoff]["hostname"]))
    nclsDynColdPlanets = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] <= cutoff]["hostname"]))
 
    nDynHotSE = len(np.unique(superEarths.loc[superEarths["st_met"] > cutoff]["hostname"]))
    nDynColdSE = len(np.unique(superEarths.loc[superEarths["st_met"] <= cutoff]["hostname"]))

    nDynHotSECompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] > cutoff]["hostname"]))
    nDynColdSECompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] <= cutoff]["hostname"]))

    nDynHotCompanions = len(np.unique(outerCompanions.loc[outerCompanions["st_met"] > cutoff]["hostname"]))
    nDynColdCompanions = len(np.unique(outerCompanions.loc[outerCompanions["st_met"] <= cutoff]["hostname"]))
    nDynHotInner = len(np.unique(planets.loc[planets["st_met"] > cutoff]["hostname"]))
    nDynColdInner = len(np.unique(planets.loc[planets["st_met"] <= cutoff]["hostname"]))

    nDynHotGG = len(np.unique(gasGiants.loc[gasGiants["st_met"] > cutoff]["hostname"]))
    nDynColdGG = len(np.unique(gasGiants.loc[gasGiants["st_met"] <= cutoff]["hostname"]))
    nDynHotGGCompanions = len(np.unique(ggCompanions.loc[ggCompanions["st_met"] > cutoff]["hostname"]))
    nDynColdGGCompanions = len(np.unique(ggCompanions.loc[ggCompanions["st_met"] <= cutoff]["hostname"]))

    occurRateCLSDynHot = stats.beta.median(nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)
    occurRateCLSDynCold = stats.beta.median(nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)
    occurRateSEDynHot = stats.beta.median(nDynHotSECompanions + 1, nDynHotSE - nDynHotSECompanions + 1)
    occurRateSEDynCold = stats.beta.median(nDynColdSECompanions + 1, nDynColdSE - nDynColdSECompanions + 1)
    occurRateCompDynHot = stats.beta.median(nDynHotCompanions + 1, nDynHotInner - nDynHotCompanions + 1)
    occurRateCompDynCold = stats.beta.median(nDynColdCompanions + 1, nDynColdInner - nDynColdCompanions +1)
    occurRateGGDynHot = stats.beta.median(nDynHotGGCompanions + 1, nDynHotGG - nDynHotGGCompanions + 1)
    occurRateGGDynCold = stats.beta.median(nDynColdGGCompanions + 1, nDynColdGG - nDynColdGGCompanions + 1)
    
    if occurRateSEDynHot > occurRateCLSDynHot:
        sigmaCLSDynHot = stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[1] - occurRateCLSDynHot
        sigmaSEDynHot = occurRateSEDynHot - stats.beta.interval(0.68, nDynHotSECompanions, nDynHotSE - nDynHotSECompanions + 1)[0]
    else:
        sigmaCLSDynHot = occurRateCLSDynHot - stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[0]
        sigmaSEDynHot = stats.beta.interval(0.68, nDynHotSECompanions, nDynHotSE - nDynHotSECompanions + 1)[1] - occurRateSEDynHot
    
    if occurRateSEDynCold > occurRateCLSDynCold:
        sigmaCLSDynCold = stats.beta.interval(0.68,nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)[1] - occurRateCLSDynCold
        sigmaSEDynCold= occurRateSEDynCold - stats.beta.interval(0.68, nDynColdSECompanions, nDynColdSE - nDynColdSECompanions + 1)[0]
    else:
        sigmaCLSDynCold = occurRateCLSDynCold -  stats.beta.interval(0.68,nclsDynColdStars + 1, nclsDynColdPlanets - nclsDynColdPlanets + 1)[0]
        sigmaSEDynCold= stats.beta.interval(0.68, nDynColdSECompanions, nDynColdSE - nDynColdSECompanions + 1)[1] - occurRateSEDynCold

    if occurRateCompDynHot > occurRateCLSDynHot:
        sigmaCLSDynHot = stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[1] - occurRateCLSDynHot
        sigmaCompDynHot = occurRateCompDynHot - stats.beta.interval(0.68, nDynHotCompanions, nDynHotInner - nDynHotCompanions + 1)[0]
    else:
        sigmaCLSDynHot = occurRateCLSDynHot - stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[0]
        sigmaCompDynHot = stats.beta.interval(0.68, nDynHotCompanions, nDynHotInner - nDynHotCompanions + 1)[1] - occurRateCompDynHot
    
    if occurRateCompDynCold > occurRateCLSDynCold:
        sigmaCLSDynCold = stats.beta.interval(0.68,nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)[1] - occurRateCLSDynCold
        sigmaCompDynCold= occurRateCompDynCold - stats.beta.interval(0.68, nDynColdCompanions, nDynColdInner - nDynColdCompanions + 1)[0]
    else:
        sigmaCLSDynCold = occurRateCLSDynCold -  stats.beta.interval(0.68,nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)[0]
        sigmaCompDynCold= stats.beta.interval(0.68, nDynColdCompanions, nDynColdInner - nDynColdCompanions + 1)[1] - occurRateCompDynCold

    if occurRateGGDynHot > occurRateCLSDynHot:
        sigmaCLSDynHot = stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[1] - occurRateCLSDynHot
        sigmaGGDynHot = occurRateGGDynHot - stats.beta.interval(0.68, nDynHotGGCompanions, nDynHotGG - nDynHotGGCompanions + 1)[0]
    else:
        sigmaCLSDynHot = occurRateCLSDynHot - stats.beta.interval(0.68,nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1)[0]
        sigmaGGDynHot = stats.beta.interval(0.68, nDynHotGGCompanions, nDynHotGG - nDynHotGGCompanions + 1)[1] - occurRateGGDynHot
    
    if occurRateGGDynCold > occurRateCLSDynCold:
        sigmaCLSDynCold = stats.beta.interval(0.68,nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)[1] - occurRateCLSDynCold
        sigmaGGDynCold= occurRateGGDynCold - stats.beta.interval(0.68, nDynColdGGCompanions, nDynColdGG - nDynColdGGCompanions + 1)[0]
    else:
        sigmaCLSDynCold = occurRateCLSDynCold -  stats.beta.interval(0.68,nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1)[0]
        sigmaGGDynCold= stats.beta.interval(0.68, nDynColdGGCompanions, nDynColdGG - nDynColdGGCompanions + 1)[1] - occurRateGGDynCold
    

    print(occurRateSEDynHot, occurRateCLSDynHot)
    
    SEDynHotEnhancement[i] = (occurRateSEDynHot - occurRateCLSDynHot)/np.sqrt(sigmaSEDynHot*sigmaSEDynHot + sigmaCLSDynHot*sigmaCLSDynHot)
    SEDynColdEnhancement[i] = (occurRateSEDynCold - occurRateCLSDynCold)/np.sqrt(sigmaSEDynCold*sigmaSEDynCold + sigmaCLSDynCold*sigmaCLSDynCold)

    CompDynHotEnhancement[i] = (occurRateCompDynHot - occurRateCLSDynHot)/np.sqrt(sigmaCompDynHot*sigmaCompDynHot + sigmaCLSDynHot*sigmaCLSDynHot)
    CompDynColdEnhancement[i] = (occurRateCompDynCold - occurRateCLSDynCold)/np.sqrt(sigmaCompDynCold*sigmaCompDynCold + sigmaCLSDynCold*sigmaCLSDynCold)

    GGDynHotEnhancement[i] = (occurRateGGDynHot - occurRateCLSDynHot)/np.sqrt(sigmaGGDynHot*sigmaGGDynHot + sigmaCLSDynHot*sigmaCLSDynHot)
    GGDynColdEnhancement[i] = (occurRateGGDynCold - occurRateCLSDynCold)/np.sqrt(sigmaGGDynCold*sigmaGGDynCold + sigmaCLSDynCold*sigmaCLSDynCold)

   

    if i in spotlightPoints:
        if i == 22:
            cutoff = 0.023
        elif i == 24:
            cutoff = 0.041
        x = np.linspace(0,0.5,1000)
        #fig1, ax1 = plt.subplots(1,1)
        #ax1.plot(x, occurrence(x, nclsDynHotPlanets + 1, nclsDynHotStars - nclsDynHotPlanets + 1), ls = "-", color = "tab:blue", label = f"P(GG|[Fe/H] > {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nclsDynColdPlanets + 1, nclsDynColdStars - nclsDynColdPlanets + 1), ls = "-", color = "tab:orange", label = f"P(GG|[Fe/H] <= {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nDynHotSECompanions + 1, nDynHotSE - nDynHotSECompanions + 1), ls = "--", color = "tab:blue", label = f"P(GG|SE,[Fe/H] > {np.round(cutoff, decimals=1)})")
        #ax1.plot(x, occurrence(x, nDynColdSECompanions + 1, nDynColdSE - nDynColdSECompanions + 1), ls = "--", color = "tab:orange", label = f"P(GG|SE,[Fe/H] <= {np.round(cutoff, decimals=1)})")
        #ax1.legend(frameon = False)
        plt.tight_layout()
        #fig1.savefig(f"./plots/occurRateMetCutoff{np.round(cutoff, decimals=2)}.png")
        print(f"Cutoff: {np.round(cutoff, decimals=5)}")
        print(f"Metal Rich Field: {nclsDynHotPlanets}/{nclsDynHotStars}, {occurRateCLSDynHot} +/- {sigmaCLSDynHot}")
        print(f"Metal Rich GG: {nDynHotGGCompanions}/{nDynHotGG}, {occurRateGGDynHot} +/- {sigmaGGDynHot}")
        print(f"Metal Poor Field: {nclsDynColdPlanets}/{nclsDynColdStars}, {occurRateCLSDynCold} +/- {sigmaCLSDynCold}")
        print(f"Metal Poor GG: {nDynColdGGCompanions}/{nDynColdGG}, {occurRateGGDynCold} +/- {sigmaGGDynCold}")
        print(f"Metal Rich Enhancement: {GGDynHotEnhancement[i]}")
        print(f"Metal Poor Enhancement: {GGDynColdEnhancement[i]}")


print(SEDynHotEnhancement)


 


fig, ax = plt.subplots(1,1)
ax.plot(eccentricityCutoffs, SEDynHotEnhancement, ls = "", marker = "o", label = "e > x, SE Inner")
ax.plot(eccentricityCutoffs, SEDynColdEnhancement, ls = "", marker = "o", label = "e <= x, SE Inner")
ax.plot(eccentricityCutoffs, GGDynHotEnhancement, ls = "", marker = "o", label = "e > x, GG Inner")
ax.plot(eccentricityCutoffs, GGDynColdEnhancement, ls = "", marker = "o", label = "e <= x, GG Inner")
ax.legend(frameon = False)
ax.set_ylabel("Enhancement ($\sigma$)")
ax.set_xlabel('Eccentricity Cutoff ')
plt.show()

#fig.savefig("./plots/metallicityEnhancement.png")

#print(eccentricityCutoffs)

    




