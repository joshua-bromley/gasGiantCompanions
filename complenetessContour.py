import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import pickle as pk
from matplotlib import rcParams

rcParams["axes.linewidth"] = 2

rcParams["ytick.right"] = True
rcParams["ytick.direction"] = "in"
rcParams["ytick.minor.visible"] = True
rcParams["ytick.major.left"] = True
rcParams["ytick.major.right"] = True
rcParams["ytick.minor.left"] = True
rcParams["ytick.minor.right"] = True
rcParams["ytick.major.size"] = 10
rcParams["ytick.minor.size"] = 5
rcParams["ytick.major.width"] = 1
rcParams["ytick.minor.width"] = 1


rcParams["xtick.top"] = True
rcParams["xtick.direction"] = "in"
rcParams["xtick.minor.visible"] = True
rcParams["xtick.major.top"] = True
rcParams["xtick.major.bottom"] = True
rcParams["xtick.minor.top"] = True
rcParams["xtick.minor.bottom"] = True
rcParams["xtick.major.size"] = 10
rcParams["xtick.minor.size"] = 5
rcParams["xtick.major.width"] = 1
rcParams["xtick.minor.width"] = 1

def averageMap(plnames):
    prob = np.zeros((len(plnames),49,49))
    for i in range(len(plnames)):
        compProb = pk.load(open(f"./completenessMaps/{plnames[i]}.p", "rb"))
        prob[i,:,:] = compProb["prob"]
    
    cmpltAvg = np.zeros((49,49))
    for j in range(49):
        for o in range(49):
            cmpltAvg[j,o] = np.sum(prob[:,j,o])
    prob = cmpltAvg/len(plnames)
    return prob

def findCountour(prob):
    smax = np.logspace(np.log10(0.3),np.log10(30),50)
    mass = np.logspace(np.log10(0.3), np.log10(30),50)
    contour = np.ones(49)*mass[0]
    for i in range(len(prob)):
        foundContour = False
        for j in range(1,len(prob[i])):
            if prob[i][j-1] <= 0.5 and prob[i][j] >= 0.5:
                contour[i] = mass[j]
                foundContour = True
        if foundContour == False:
            if prob[i][-1] < 0.5:
                contour[i] = mass[-1]
            elif prob[i][0] > 0.5:
                contour[i] = mass[0]
    
    return contour

planets = pd.read_csv("./data/gasGiantDataComplete.csv")
planets = planets.loc[(pd.isna(planets["pl_bmassj"]) == False) & (pd.isna(planets["pl_orbsmax"]) == False) & (planets["pl_orbeccen"] > 0)]
coldJupiters = planets.loc[(planets["pl_type"] == "CJ")]
warmJupiters = planets.loc[(planets["pl_type"] == "WJ")]
hotJupiters = planets.loc[planets["pl_type"] == 'HJ']
hotSaturns = planets.loc[(planets["pl_type"] == "HS")]
subSaturns = planets.loc[(planets["pl_type"] == "CS")]
superEarths = planets.loc[planets["pl_type"] == "SE"]
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
hsCompanions = planets.loc[(planets["companion_type"] % 3 == 0)]
csCompanions = planets.loc[planets["companion_type"] % 5 == 0]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
wjCompanions = planets.loc[planets["companion_type"] % 11 == 0]
cjCompanions = planets.loc[ (planets["companion_type"] % 13 == 0)] 

#planetsList = [coldJupiters, warmJupiters, hotJupiters, seCompanions, hsCompanions, csCompanions, hjCompanions, wjCompanions, cjCompanions]
#labels = ["Cold Jupiters","Warm Jupiters", "Hot Jupiters", "SE Compnaions", "HS Companions","CS Companions", "HJ Companions", "WJ Companions", "CJ Companions"]

planetsList = [superEarths, hotSaturns, subSaturns, hotJupiters, warmJupiters, coldJupiters]
labels = ["Super Earths", "Hot Sub-Saturns", "Cold Sub-Saturns", "Hot Jupiters", "Warm Jupiters", "Cold Jupiters"]
colours = ["tab:orange", "tab:green", "tab:cyan", "tab:red", "tab:pink", "tab:purple"]


fig, ax =  plt.subplots(1,1)

for i in range(len(planetsList)):
    smax = np.logspace(np.log10(0.3),np.log10(30),49)
    cmpltMap = averageMap(np.unique(planetsList[i]["hostname"]))
    contour = findCountour(cmpltMap)
    ax.plot(smax, contour, label = labels[i], color = colours[i])

ax.legend(frameon = False, fontsize = 16)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.3,30)
ax.set_ylim(0.3,30)
ax.set_xlabel("Semi-Major Axis (AU)", fontsize = 16)
ax.set_ylabel("Mass ($M_J$)", fontsize = 16)

tickLabelSize = 12
ax.tick_params(axis = 'x', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'x', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "major", direction = "in", labelsize = tickLabelSize, pad = 10)
ax.tick_params(axis = 'y', bottom = True, top = True, which = "minor", direction = "in", labelsize = tickLabelSize, pad = 10)
plt.tight_layout()
fig.savefig("./plots/completenessContours.png")
plt.show()

        

    
