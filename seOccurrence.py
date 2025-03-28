import numpy as np
import pandas as pd
from scipy import special as sp

clsStars = pd.read_csv("./data/clsStars.csv")
clsPlanets = pd.read_csv("./data/clsPlanets2.csv", skiprows = 101)
clsPlanets = clsPlanets.loc[(clsPlanets["pl_bmasse"] > 0.5*317) & (clsPlanets["pl_orbsmax"] > 1)]

planets = pd.read_csv("./data/gasGiantData.csv")

superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[(planets["companion_type"] % 2 == 0) & (planets["pl_orbsmax"] > 1)]

nGasGiant = len(np.unique(clsPlanets["hostname"]))
nGasGiantMR = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] > 0]["hostname"]))
nGasGiantMP = len(np.unique(clsPlanets.loc[clsPlanets["st_met"] <= 0]["hostname"]))

nclsStars = len(clsStars)
nclsStarsMR = len(clsStars.loc[clsStars["[Fe/H]"] > 0])
nclsStarsMP = len(clsStars.loc[clsStars["[Fe/H]"] <= 0])

nSECompanions = len(np.unique(seCompanions["hostname"]))
nSECompanionsMR = len(np.unique(seCompanions.loc[seCompanions["st_met"] > 0]["hostname"]))
nSECompanionsMP = len(np.unique(seCompanions.loc[seCompanions["st_met"] <= 0]["hostname"]))

nSuperEarth = len(np.unique(superEarths["hostname"]))
nSuperEarthMR = len(np.unique(superEarths.loc[superEarths["st_met"] > 0]["hostname"]))
nSuperEarthMP = len(np.unique(superEarths.loc[superEarths["st_met"] <= 0]["hostname"]))

print("P(GG   | ------------ | |[Fe/H] > 0)           | |[Fe/H]<= 0)")
print("--------------------------------------------------")
print(f"----- | {nGasGiant/nclsStars} | {nGasGiantMR/nclsStarsMR} | {nGasGiantMP/nclsStarsMP}")
print(f" |SE) | {nSECompanions/nSuperEarth} | {nSECompanionsMR/nSuperEarthMR} | {nSECompanionsMP/nSuperEarthMP}")\


print(nclsStars, nclsStarsMR, nclsStarsMP)
