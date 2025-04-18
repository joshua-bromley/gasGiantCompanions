import pandas as pd
import numpy as np

systemData = pd.read_csv("./RVData/rvData.txt")
systems = systemData["hostname"].values

for sys in systems:
    try:
        if sys == "XO2S" or sys == "WASP53":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter="\t", usecols = (0,1,2)).T
        elif sys == "HD72659" or sys == "HD68988" or sys == "HD38529" or sys == "HD28185" or sys == "HD24040" or sys == "HD183263" or sys == "HD169830" or sys == "HD168443" or sys == "HD156279" or sys == "HD12661" or sys == "14Her" or sys == "HD87883" or sys == "HD8574" or sys == "HD80606" or sys == "HD50554" or sys == "HD213472" or sys == "HD209458":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (5,6,7)).T
            time = time - 2450000
        else:
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (0,1,2)).T
    except ValueError:
        print("Bad Data for system " + sys)
    