import pandas as pd
import numpy as np

systemData = pd.read_csv("./RVData/rvData.txt")

defaultSource = systemData.loc[pd.isna(systemData["data_ref"]) == True]
clsSource = systemData.loc[systemData["data_ref"] == "Rosenthal+21"]
bonomoSource = systemData.loc[(systemData["prop_ref"] == "Bonomo+17") | (systemData["data_ref"] == "Bonomo+17")]
otherSource = systemData[(pd.isna(systemData["data_ref"]) == False) & (systemData["data_ref"] != "Rosenthal+21") & (systemData["prop_ref"] != "Bonomo+17") & (systemData["data_ref"] != "Bonomo+17")]

np.savetxt("./namelists/nameListDefaultSource.txt", defaultSource["hostname"].values, fmt="%s")
np.savetxt("./namelists/namelistCLSSource.txt", clsSource["hostname"].values, fmt="%s")
np.savetxt("./namelists/namelistBonomoSource.txt", bonomoSource["hostname"].values, fmt="%s")
np.savetxt("./namelists/namelistOtherSource.txt", otherSource[["hostname", "data_ref"]].values, fmt="%s")

