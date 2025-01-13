import numpy as np
import pandas as pd
mJ = 317.8

planets = pd.read_csv("./data/gasGiantDataRaw.csv", skiprows=101)

columns = ["pl_name", "hostname", "sy_snum", "sy_pnum", "pl_orbper", "pl_orbsmax", "pl_bmasse", "pl_bmassj", "pl_orbeccen","st_spectype", "st_teff", "st_rad", "st_mass", "st_met"]
planetTypes = np.empty_like(planets["pl_name"])
ggCompanion = np.empty_like(planets["pl_name"])
ssCompanion = np.empty_like(planets["pl_name"])
seCompanion = np.empty_like(planets["pl_name"])

for i in range(len(planets)):
    if pd.isnull(planets.iloc[i]["pl_orbsmax"]):
        smax = (planets.iloc[i]["st_mass"]/(planets.iloc[i]["pl_orbper"]/365)**2)**(1/3)
        planets.at[i, "pl_orbsmax"] = smax
        
    if planets.iloc[i]["pl_bmasse"] < 20:
        planetTypes[i] = "SE"
    elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 0.1:
        planetTypes[i] = "HS"
    elif planets.iloc[i]["pl_bmasse"] < 0.5*mJ and planets.iloc[i]["pl_orbsmax"] > 0.1:
        planetTypes[i] = "CS"
    elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 0.1:
        planetTypes[i] = "HJ"
    elif planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["pl_orbsmax"] < 1:
        planetTypes[i] = "WJ"
    else:
        planetTypes[i] = "CJ"

    ggCompanion[i] = 0
    ssCompanion[i] = 0
    seCompanion[i] = 0
    if planets.iloc[i]["pl_bmasse"] >= 0.5*mJ and planets.iloc[i]["sy_pnum"] > 1:
        hostname = planets.iloc[i]["hostname"]
        companions = planets.loc[planets["hostname"] == hostname]
        for j in range(len(companions)):
            if companions.iloc[j]["pl_orbsmax"] < planets.iloc[i]["pl_orbsmax"]:
                if companions.iloc[j]["pl_bmasse"] < 20:
                    seCompanion[i] = 1
                elif companions.iloc[j]["pl_bmasse"] < 0.5*mJ:
                    ssCompanion[i] = 1
                else:
                    ggCompanion[i] = 1
        

planets = planets[columns]
planets.insert(4, "pl_type", planetTypes)
planets.insert(15, "gg_companion", ggCompanion)
planets.insert(15, "ss_companion", ssCompanion)
planets.insert(15, "se_companion", seCompanion)

planets.to_csv("./data/gasGiantData.csv")
    
