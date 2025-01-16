import numpy as np
import pandas as pd
mJ = 317.8

planets = pd.read_csv("./data/gasGiantDataRaw.csv", skiprows=101)

columns = ["pl_name", "hostname","pl_type", "sy_snum", "sy_pnum", "pl_orbper", "pl_orbsmax","pl_rade","pl_radj", "pl_bmasse", "pl_bmassj", "pl_orbeccen","st_spectype", "st_teff", "st_rad", "st_mass", "st_met"]
planetTypes = np.empty_like(planets["pl_name"])
companionType = np.ones_like(planets["pl_bmasse"])
suspiciousData = np.zeros_like(planets["pl_bmasse"])

for i in range(len(planets)):
    if pd.isnull(planets.iloc[i]["pl_orbsmax"]):
        smax = (planets.iloc[i]["st_mass"]/(planets.iloc[i]["pl_orbper"]/365)**2)**(1/3)
        planets.at[i, "pl_orbsmax"] = smax

    limitFlag = planets.iloc[i]["pl_orbperlim"] + planets.iloc[i]["pl_orbsmaxlim"] + planets.iloc[i]["pl_bmasselim"] + planets.iloc[i]["pl_orbeccenlim"]
    if pd.isnull(planets.iloc[i]["pl_bmasse"]) or pd.isnull(planets.iloc[i]["pl_orbsmax"]) or limitFlag > 0:
        suspiciousData[i] = 1

    if not pd.isnull(planets.iloc[i]["pl_bmasse"]):  
        if not pd.isnull(planets.iloc[i]["pl_orbsmax"]):  
            if planets.iloc[i]["pl_bmasse"] < 20:
                planetTypes[i] = "SE"
            elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 0.1:
                planetTypes[i] = "HS"
            elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbsmax"] >= 0.1:
                planetTypes[i] = "CS"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 0.1:
                planetTypes[i] = "HJ"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 1:
                planetTypes[i] = "WJ"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbsmax"] >= 1:
                planetTypes[i] = "CJ"
        elif not pd.isnull(planets.iloc[i]["pl_orbper"]):
            if planets.iloc[i]["pl_bmasse"] < 20:
                planetTypes[i] = "SE"
            elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbper"] < 10:
                planetTypes[i] = "HS"
            elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbper"] >= 10:
                planetTypes[i] = "CS"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbper"] < 10:
                planetTypes[i] = "HJ"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbper"] < 300:
                planetTypes[i] = "WJ"
            elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbper"] >= 300:
                planetTypes[i] = "CJ"
        else:
            planetTypes[i] = "NA"
    elif not pd.isnull(planets.iloc[i]["pl_rade"]):
        if not pd.isnull(planets.iloc[i]["pl_orbsmax"]):  
            if planets.iloc[i]["pl_rade"] < 4:
                planetTypes[i] = "SE"
            elif planets.iloc[i]["pl_rade"] < 9 and planets.iloc[i]["pl_orbsmax"] < 0.1:
                planetTypes[i] = "HS"
            elif planets.iloc[i]["pl_rade"] < 9 and planets.iloc[i]["pl_orbsmax"] >= 0.1:
                planetTypes[i] = "CS"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbsmax"] < 0.1:
                planetTypes[i] = "HJ"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbsmax"] < 1:
                planetTypes[i] = "WJ"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbsmax"] >= 1:
                planetTypes[i] = "CJ"
        elif not pd.isnull(planets.iloc[i]["pl_orbper"]):
            if planets.iloc[i]["pl_rade"] < 4:
                planetTypes[i] = "SE"
            elif planets.iloc[i]["pl_rade"] < 9 and planets.iloc[i]["pl_orbper"] < 10:
                planetTypes[i] = "HS"
            elif planets.iloc[i]["pl_rade"] < 9 and planets.iloc[i]["pl_orbper"] >= 10:
                planetTypes[i] = "CS"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbper"] < 10:
                planetTypes[i] = "HJ"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbper"] < 300:
                planetTypes[i] = "WJ"
            elif planets.iloc[i]["pl_rade"] >= 9 and planets.iloc[i]["pl_orbper"] >= 300:
                planetTypes[i] = "CJ"
        else:
            planetTypes[i] = "NA"
    else:
        planetTypes[i] = "NA"

planets.insert(4, "pl_type", planetTypes)

for i in range(len(planets)):
    if planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["sy_pnum"] > 1:
        hostname = planets.iloc[i]["hostname"]
        companions = planets.loc[planets["hostname"] == hostname]
        for j in range(len(companions)):
            if (companions.iloc[j]["pl_orbsmax"] < planets.iloc[i]["pl_orbsmax"]) or (companions.iloc[j]["pl_orbper"] < planets.iloc[i]["pl_orbper"]):
                if companions.iloc[j]["pl_type"] == "SE":
                    companionType[i] *= 2
                elif companions.iloc[j]["pl_type"] == "HS":
                    companionType[i] *= 3
                elif companions.iloc[j]["pl_type"] == "CS":
                    companionType[i] *= 5
                elif companions.iloc[j]["pl_type"] == "HJ":
                    companionType[i] *= 7
                elif companions.iloc[j]["pl_type"] == "WJ":
                    companionType[i] *= 11
                elif companions.iloc[j]["pl_type"] == "CJ":
                    companionType[i] *= 13        

"""
Companion Type is a product of primes so that I can identify the categories a companion belongs to with a single number
2: Super Earths
3: Hot Saturns
5: Cold Saturns
7: Hot Jupiters
11: Warm Jupiters
13: Cold Jupiters
"""

planets = planets[columns]
planets.insert(17, "companion_type", companionType)
planets.insert(18, "pl_sus_data", suspiciousData)


planets.to_csv("./data/gasGiantData.csv")

print(np.sum(suspiciousData)/len(suspiciousData))
    
