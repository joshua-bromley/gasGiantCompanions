import pandas as pd
import numpy as np

planets = pd.read_csv("./data/gasGiantComplete2.csv")
metallicities = pd.read_csv("./data/metallicities.csv")

for i in range(len(planets)):
    if pd.isna(planets.iloc[i]["st_met"]):
        hostname = planets.iloc[i]["hostname"]
        metallicity = metallicities.loc[metallicities["hostname"] == hostname]
        if len(metallicity) > 0:
            met = metallicity.iloc[0]["st_met"]
            planets.at[i, "st_met"] = met

planets.to_csv("./data/gasGiantDataComplete2.csv")