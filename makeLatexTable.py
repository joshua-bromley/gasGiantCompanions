import numpy as np
import pandas as pd

planets = pd.read_csv("./data/gasGiantDataCompleteTable.csv")
planets = planets.loc[(planets["st_mass"] > 0.6) & (planets["st_mass"] < 1.6)]
references = np.unique(planets["pl_refname"])
hostnames = np.unique(planets["hostname"])
rows = []

for i in range(len(hostnames)):
    hostname = hostnames[i]
    system = planets.loc[planets["hostname"] == hostname].sort_values(by="pl_bmassj", ascending = False)
    nPlanets = "{:.0f}".format(len(system))
    largestPlanet = system.iloc[0]
    stMass = "{:.2f}".format(largestPlanet["st_mass"])
    stMet = "{:.2f}".format(largestPlanet["st_met"])
    baseline = "{:.0f}".format(largestPlanet["baseline"])
    rvDataPoints = "{:.0f}".format(largestPlanet["num_obs"])
    reference = "("+"{:.0f}".format(np.argwhere(references == largestPlanet["pl_refname"])[0][0]+1)+")"
    planetTypes = np.zeros(4)
    planetTypeStr = "" 
    for j in range(len(system)):
        plType = system.iloc[j]["pl_type"]
        match plType:
            case "CJ":
                if planetTypes[0] == 0:
                    planetTypeStr += "CJ, "
                planetTypes[0] += 1
            case "WJ":
                if planetTypes[0] == 0:
                    planetTypeStr += "WJ, "
                planetTypes[0] += 1
            case "HJ":
                if planetTypes[1] == 0:
                    planetTypeStr += "HJ, "
                planetTypes[1] += 1
            case "CS":
                if planetTypes[2] == 0:
                    planetTypeStr += "CS, "
                planetTypes[2] += 1
            case "HS":
                if planetTypes[2] == 0:
                    planetTypeStr += "HS, "
                planetTypes[2] += 1
            case "SE":
                if planetTypes[3] == 0:
                    planetTypeStr += "SE, "
                planetTypes[3] += 1
    planetTypeStr = planetTypeStr[:-2]
    hasGasGiant = "N"
    if len(system) > 1 and (largestPlanet["pl_type"] == "CJ" or largestPlanet["pl_type"] == "WJ"):
        hasGasGiant = "Y"
    tableEntry = hostname + " & " + nPlanets + " & " + hasGasGiant + " & " + planetTypeStr + " & " + stMass + " & " + stMet + " & " + baseline + " & " + rvDataPoints + " & " + reference + " \\\\"
    rows.append(tableEntry)

np.savetxt("./latexTable.txt", rows, fmt = "%s")
                



    

