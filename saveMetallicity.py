import numpy as np
import pandas as pd

planets = pd.read_csv("./data/gasGiantData.csv")
metallicities = planets[["pl_name", "st_met"]]
metallicities.to_csv("./data/metallicities.csv")

