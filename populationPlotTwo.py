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

planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]

hosts = np.unique(planets["hostname"])

fig, ax = plt.subplots(1,1, figsize = (6,4))

for host in hosts:
    system = planets.loc[planets["hostname"] == host].sort_values(by = "pl_orbsmax")
    ax.plot(system["pl_orbsmax"].values, system["pl_bmassj"].values, color = "tab:blue", alpha = 0.3)

for i in range(len(planets)):
    mass = planets.iloc[i]["pl_bmassj"]
    smax = planets.iloc[i]["pl_orbsmax"]
    colour = "tab:purple"
    if planets.iloc[i]["pl_type"] == "HJ":
        colour = "tab:red"
    if planets.iloc[i]["pl_type"] == "SE":
        colour = "tab:orange"
    if planets.iloc[i]["pl_type"] == "HS" or planets.iloc[i]["pl_type"] == "CS":
        colour = "tab:green"
    ax.plot(smax, mass, color = colour, marker = "o")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Semi Major Axis (AU)", fontsize = 16)
ax.set_ylabel("Mass ($M_J$)", fontsize = 16)
tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()

fig.savefig("./plots/cmpltCorrPopulation.pdf")
plt.show()

