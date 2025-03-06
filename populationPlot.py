import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import special as sp

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

planets  = pd.read_csv("./data/gasGiantData.csv")
hosts = np.unique(planets["hostname"])

ssMass = [0.815/318, 1/318, 318/318, 95/318, 14.5/318, 17/318]
ssSmax = [ 0.723, 1, 5.203, 9.537, 19.191, 30.069]
ssPer = [ 224.701, 365.26, 4332.59, 10759.22, 30688.5, 60182]

fig, ax = plt.subplots(1,1, figsize = (6,4))
for hostname in hosts:
    system = planets.loc[planets["hostname"] == hostname].sort_values("pl_orbsmax")
    mass = system["pl_bmassj"].values
    smax = system["pl_orbsmax"].values
    period = system["pl_orbper"].values
    ax.plot(np.log10(period), np.log10(mass), marker = ".", ls = " ", color = "tab:blue")
    #ax.plot(period, mass, marker = ".", color = "tab:blue", alpha = 0.3)

ax.plot(np.log10(ssPer), np.log10(ssMass), marker = ".", ls = "", color = "tab:pink")

xTicks = [-1,0,1,2,3,4,5]
xMinorTicks = [-0.5,0.5,1.5,2.5,3.5,4.5]
yTicks = [-3,-2,-1,0,1,2]
yMinorTicks = [-2.5,-1.5,-0.5,0.5,1.5]

xTickLabels = ["$10^{-1}$", "1", "10", "$10^2$", "$10^3$", "$10^4$", "$10^5$"]
yTickLabels = ["$10^{-3}$","$10^{-2}$","$10^{-1}$", "1", "10", "$10^2$"]

ax.set_xticks(xTicks)
ax.set_xticks(xMinorTicks, minor=True)
ax.set_yticks(yTicks)
ax.set_yticks(yMinorTicks, minor=True)
ax.set_xticklabels(xTickLabels)
ax.set_yticklabels(yTickLabels)
ax.set_xlim(-1,5)
ax.set_ylim(-3,2)

ax.set_xlabel("Period (Days)", fontsize = 16)
ax.set_ylabel("Mass (Jupiter Masses)", fontsize = 16)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()

fig.savefig("./plots/planetSample.png")
fig.savefig("./plots/planetSample.pdf") 

gasGiants = planets.loc[(planets["pl_type"] == "HJ") | (planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
seCompanions = planets.loc[planets["companion_type"] % 2 == 0]
wcjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]
hjCompanions = planets.loc[(planets["companion_type"] % 7 == 0)]

fig2, ax2 = plt.subplots(1,1)
ax2.ecdf(gasGiants["pl_orbeccen"].dropna(), label = "All Gas Giants")
ax2.ecdf(seCompanions["pl_orbeccen"].dropna(), label = "Super Earth Companions")
ax2.ecdf(wcjCompanions["pl_orbeccen"].dropna(), label = "Warm/Cold Jupiter Companions")
ax2.ecdf(hjCompanions["pl_orbeccen"].dropna(), label = "Hot Jupiter Companions")
ax2.legend(frameon = False)

ax2.set_xlabel("Eccentricity", fontsize = 16)
ax2.set_ylabel("Proportion (CDF)", fontsize = 16)

tickLabelSize = 12
ax2.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax2.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax2.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax2.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)

ax2.set_xlim(-0.05, 1.05)

plt.tight_layout()
fig2.savefig("./plots/eccenDist.png")
fig2.savefig("./plots/eccenDist.pdf")

superEarths = planets.loc[planets["pl_type"] == "SE"]
warmColdJupiters = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == 'CJ')]

nCLSMetalRichStars = 395
nCLSMetalPoorStars = 324

nCLSMetalRichGG = 49
nCLSMetalPoorGG = 13

nMRseCompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] >= 0]["hostname"]))
nMPseCompanions = len(np.unique(seCompanions.loc[seCompanions["st_met"] < 0]["hostname"]))
nMRsuperEarths = len(np.unique(superEarths.loc[superEarths["st_met"] >= 0]["hostname"]))
nMPsuperEarths = len(np.unique(superEarths.loc[superEarths["st_met"] < 0]["hostname"]))

nMRwcjCompanions = len(np.unique(wcjCompanions.loc[wcjCompanions["st_met"] >= 0]["hostname"]))
nMPwcjCompanions = len(np.unique(wcjCompanions.loc[wcjCompanions["st_met"] < 0]["hostname"]))
nMRwarmColdJupiters = len(np.unique(warmColdJupiters.loc[warmColdJupiters["st_met"] >= 0]["hostname"]))
nMPwarmColdJupiters = len(np.unique(warmColdJupiters.loc[warmColdJupiters["st_met"] < 0]["hostname"]))



def occurance(x, a, b):
    return (1/sp.beta(a,b))* x**a * (1-x)**b

x = np.linspace(0,0.8,1000)
fig3, ax3 = plt.subplots(1,1, figsize = (6,4))
ax3.plot(x, occurance(x, nCLSMetalRichGG, nCLSMetalRichStars - nCLSMetalRichGG), color = "tab:blue", label = "P(GG|[Fe/H] >= 0)")
ax3.plot(x, occurance(x, nCLSMetalPoorGG, nCLSMetalPoorStars - nCLSMetalPoorGG), color = "tab:blue", ls = "--", label = "P(GG|[Fe/H] < 0)")
ax3.plot(x, occurance(x, nMRseCompanions, nMRsuperEarths - nMRseCompanions), color = "tab:orange", label = "P(GG|SE,[Fe/H] >= 0)")
ax3.plot(x, occurance(x, nMPseCompanions, nMPsuperEarths - nMPseCompanions), color ="tab:orange", ls = "--", label = "P(GG|SE, [Fe/H] < 0)")
ax3.plot(x, occurance(x, nMRwcjCompanions, nMRwarmColdJupiters - nMRwcjCompanions), color = "tab:green", label = "P(GG|GG,[Fe/H] >= 0)")
ax3.plot(x, occurance(x, nMPwcjCompanions, nMPwarmColdJupiters - nMPwcjCompanions), color ="tab:green", ls = "--", label = "P(GG|GG, [Fe/H] < 0)")
ax3.legend(frameon = False, fontsize = tickLabelSize)

ax3.set_xlabel("P(GG)", fontsize = 16)
ax3.set_ylabel("PDF", fontsize = 16)

ax3.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax3.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax3.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax3.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)

ax3.set_xlim(-0.02,0.82)
ax3.set_ylim(-0.1,2.7)
plt.tight_layout()
fig3.savefig("./plots/metProbs.png")
fig3.savefig("./plots/metProbs.pdf")
 

