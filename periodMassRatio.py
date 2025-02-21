import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

planets = pd.read_csv("./data/gasGiantData.csv").drop_duplicates(subset = ["pl_name"])

planetPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])

hostnames = np.unique(planets["hostname"])

for i in range(len(hostnames)):
    hostname = hostnames[i]
    system = planets.loc[planets["hostname"] == hostname].sort_values(by="pl_orbper")
    for j in range(1,len(system)):
        planetA = system.iloc[j]
        planetB = system.iloc[j-1]
        planetAName = planetA["pl_name"]
        planetBName = planetB["pl_name"]
        planetAMass = planetA["pl_bmasse"]
        planetAPer = planetA["pl_orbper"]
        planetAeccen = planetA["pl_orbeccen"]
        planetBMass = planetB["pl_bmasse"]
        planetBPer = planetB["pl_orbper"]
        planetBeccen = planetB["pl_orbeccen"]
        stellarMass = planetA["st_mass"]
        stellarMet = planetA["st_met"]
        massRatio = planetAMass/planetBMass
        perRatio = planetAPer/planetBPer
        planetPairs.loc[len(planetPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]


highMassPairs = planetPairs.loc[planetPairs["pla_mass"] > 0.5*317]

periodRatios = planetPairs["per_ratio"].values
massRatios = planetPairs["mass_ratio"].values
outerEccen = planetPairs["pla_eccen"].values

fig, ax = plt.subplots(1,1)
ax.plot(periodRatios, massRatios, ls = "", marker = "o")
ax.set_xscale("log")
ax.set_yscale("log")


plt.show()

        
