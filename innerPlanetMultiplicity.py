import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import special as sp

planets = pd.read_csv("./data/gasGiantData.csv")

smallPlanets = planets.loc[(planets["pl_type"] == "SE") &(planets["st_met"] > 0)]# | (planets["pl_type"] == "HS") | (planets["pl_type"] == "CS")]
spCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["st_met"] > 0)]# | (planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]

multiSmallPlanets = []
mspCompanions = []
for i in range(len(smallPlanets["pl_name"])):
    hostname = smallPlanets.iloc[i]["hostname"]
    smallCompanions = smallPlanets.loc[smallPlanets["hostname"] == hostname]
    if len(smallCompanions["pl_name"]) > 1:
        multiSmallPlanets.append(hostname)
    if hostname in np.array(spCompanions["hostname"]):
        mspCompanions.append(hostname)

multiSmallPlanets = np.unique(multiSmallPlanets)
mspCompanions = np.unique(mspCompanions)

print(len(multiSmallPlanets))
print(len(mspCompanions))

def occurrance(x,a,b):
    return (1/sp.beta(a,b)) * x**a * (1-x)**(b)

nSmallPlanets = len(np.unique(smallPlanets["hostname"]))
nSPCompanions = len(np.unique(spCompanions["hostname"]))

x = np.linspace(0,1,1000)
fig, ax = plt.subplots(1,1)
probGGSP = occurrance(x, nSPCompanions, nSmallPlanets - nSPCompanions)
probGGMSP = occurrance(x, len(mspCompanions), len(multiSmallPlanets) - len(mspCompanions))
ax.plot(x, probGGSP, label = "All Inner Planet Companions")
ax.plot(x, probGGMSP, label = "Multiple Inner Companions")
ax.legend(frameon = False)
fig.savefig("./plots/innerPlanetMultiplicity.png")