import pandas as pd
import numpy as np
import pickle as pk

systemData = pd.read_csv("./RVData/rvData.txt")
systems = systemData["hostname"].values
systems = ["HD6061", "HD63433", "HD63935", "HD77946", "HD93963", "HIP8152", "HIP9618", "HIP97166", "TOI1180", "TOI1184", "TOI1194", "TOI1249", "TOI1272", "TOI1279", "TOI1288", "TOI1294", "TOI1296", "TOI1298", "TOI1386", "TOI1410", "TOI1439", "TOI1443", "TOI1444", "TOI1451", "TOI1472", "TOI1601", "TOI1691", "TOI1710", "TOI1723", "TOI1736", "TOI1742", "TOI1751", "TOI1794", "TOI1798", "TOI1799", "TOI1807", "TOI1823", "TOI1824", "TOI1842", "TOI1898", "TOI2019", "TOI2076", "TOI2128", "TOI329", "TOI480", "TOI669", "WASP156"]


for sys in systems:
    try:
        if sys == "XO2S" or sys == "WASP53":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter="\t", usecols = (0,1,2)).T
        elif sys == "HD72659" or sys == "HD68988" or sys == "HD38529" or sys == "HD28185" or sys == "HD24040" or sys == "HD183263" or sys == "HD169830" or sys == "HD168443" or sys == "HD156279" or sys == "HD12661" or sys == "14Her" or sys == "HD87883" or sys == "HD8574" or sys == "HD80606" or sys == "HD50554" or sys == "HD213472" or sys == "HD209458" or sys == "HD189733" or sys == "HD136925" or sys == "HD3765" or sys == "HD3651" or sys == "HD164922" or sys == "HD156668" or sys == "HD107148" or sys == "16Cyg" or sys == "HD104067" or sys == "HD10697" or sys =="HD114729" or sys == "HD117207" or sys == "HD126614" or sys == "HD13931" or sys == "HD154345" or sys == "HD16141" or sys == "HD170469" or sys == "HD178911B" or sys == "HD181234" or sys == "HD192263" or sys == "HD195019" or sys == "HD210277" or sys == "HD217107" or sys == "HD218566":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (5,6,7)).T
            time = time - 2450000
        elif sys == "HD222582" or sys == "HD26161" or sys == "HD34445" or sys == "HD40979" or sys == "HD4208" or sys == "HD45350" or sys == "HD66428" or sys == "HD82943" or sys == "HD99109" or sys == "HD168746" or sys == "HD49674":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (5,6,7)).T
            time = time - 2450000
        elif sys == "TOI238":
             time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (0,1,8)).T
        elif sys == "HD6061" or sys == "HD63433" or sys == "HD63935" or sys == "HD77946" or sys == "HD93963" or sys == "HIP8152" or sys == "HIP9618" or sys == "HIP97166" or sys == "TOI1180" or sys == "TOI1184" or sys == "TOI1194" or sys == "TOI1249" or sys == "TOI1296" or sys == "TOI1272" or sys == "TOI1279" or sys == "TOI1288" or sys == "TOI1294" or sys == "TOI1298" or sys == "TOI1386" or sys == "TOI1410" or sys == "TOI1439" or sys == "TOI1443" or sys == "TOI1444" or sys == "TOI1451" or sys == "TOI1472" or sys == "TOI1601" or sys == "TOI1691":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (2,3,4)).T
            time = time - 2450000
        elif sys == "TOI1710" or sys == "TOI1723" or sys == "TOI1736" or sys == "TOI1742" or sys == "TOI1751" or sys == "TOI1794" or sys == "TOI1798" or sys == "TOI1799" or sys == "TOI1807" or sys == "TOI1823" or sys == "TOI1824" or sys == "TOI1842" or sys == "TOI1898" or sys == "TOI2019" or sys == "TOI2076" or sys == "TOI2128" or sys == "TOI329" or sys == "TOI480" or sys == "TOI669" or sys == "WASP156":
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (2,3,4)).T
            time = time - 2450000
        else:
            time, rv, err = np.loadtxt(f"./RVData/rv{sys}.txt", delimiter=",", usecols = (0,1,2)).T
    except ValueError as e:
        print("Bad Data for system " + sys)
        print(e)
    except FileNotFoundError as e:
        print("No RV data for system " + sys)
    try:
        compProb = pk.load(open(f"./completenessMaps/{sys}.p", "rb"))
    except FileNotFoundError:
        print("No Completeness Map for " + sys)
    