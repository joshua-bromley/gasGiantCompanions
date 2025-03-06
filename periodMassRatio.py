import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

planets = pd.read_csv("./data/gasGiantData.csv").drop_duplicates(subset = ["pl_name"])

planetPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
hjPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])
wcjPairs = pd.DataFrame(columns=["hostname","planet_a_name", "planet_b_name", "pla_per", "plb_per","per_ratio", "pla_mass", "plb_mass", "mass_ratio", "pla_eccen", "plb_eccen","st_mass","st_met"])


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
        if planetB["pl_type"] == "HJ" or planetA["pl_type"] == "HJ":
            hjPairs.loc[len(hjPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
        elif (planetB["pl_type"] == "WJ" or planetB["pl_type"] == "CJ") or (planetA["pl_type"] == "WJ" or planetA["pl_type"] == "CJ"):
            wcjPairs.loc[len(wcjPairs)] = [hostname, planetAName, planetBName, planetAPer, planetBPer, perRatio, planetAMass, planetBMass, massRatio, planetAeccen, planetBeccen, stellarMass, stellarMet]
 

highMassPairs = planetPairs.loc[planetPairs["pla_mass"] > 0.5*317]

periodRatios = planetPairs["per_ratio"].values
massRatios = planetPairs["mass_ratio"].values
outerEccen = planetPairs["pla_eccen"].values

hjPerRatios = hjPairs["per_ratio"].values
hjMassRatios = hjPairs["mass_ratio"].values

wcjPerRatios = wcjPairs["per_ratio"].values
wcjMassRatios = wcjPairs["mass_ratio"].values

fig, ax = plt.subplots(1,1)
ax.plot(hjPairs["plb_mass"].values, hjPairs["pla_mass"].values, ls = "", marker = "o")
ax.plot(wcjPairs["plb_mass"].values, wcjPairs["pla_mass"].values, ls = "", marker = "o")
for i in range(len(hjPairs)):
    ax.text(hjPairs["plb_mass"].values[i], hjPairs["pla_mass"].values[i], hjPairs.iloc[i]["hostname"])
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Inner Planet Mass (M_E)")
ax.set_ylabel('Outer Planet Mass (M_E)')
#fig.savefig("./plots/perMassRatio.png")



print(np.mean(hjMassRatios))
print(np.mean(wcjMassRatios))

        

plt.show()