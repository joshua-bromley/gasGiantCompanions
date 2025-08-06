import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams

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

planets = pd.read_csv("./data/gasGiantComplete2.csv")
companions = planets.loc[planets["companion_type"] > 1]

planets = planets.drop_duplicates(subset = "hostname")
companions = companions.drop_duplicates(subset = "hostname")

fig, ax, = plt.subplots(1,2, figsize = (8,4))
massBins = np.arange(0,3.1,0.25)
metBins = np.arange(-1,0.61,0.1)

ax[0].hist(planets["st_mass"], bins = massBins, label = "All Systems")
ax[0].hist(companions["st_mass"], bins = massBins, alpha = 0.8, label = "Outer Gas Giants")
ax[1].hist(planets["st_met"], bins = metBins, label = "All Systems")
ax[1].hist(companions["st_met"], bins = metBins, label = "Outer Gas Giants", alpha = 0.8)

ax[0].set_ylabel("Count", fontsize = 16)
ax[0].set_xlabel("Stellar Mass ($M_\odot$)", fontsize = 16)
ax[1].set_xlabel("Metallicity ([Fe/H])", fontsize = 16)
ax[1].legend(frameon = False, fontsize = 12)

xTicksLeft = np.arange(0,3.1,0.5)
yTicksLeft = np.arange(0,251,50)
ax[0].set_xticks(xTicksLeft)
ax[0].set_yticks(yTicksLeft)
ax[0].set_xlim(0,3.0)
ax[0].set_ylim(0,250)

xTicksRight = np.arange(-1,0.51,0.5)
yTicksRight = np.arange(0,121,20)
ax[1].set_xticks(xTicksRight)
ax[1].set_yticks(yTicksRight)
ax[1].set_xlim(-1.1,0.6)
ax[1].set_ylim(0,120)

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

fig.savefig("./plots/stellarPropertiesAll.png")
plt.show()
