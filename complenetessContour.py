import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import pickle as pk

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
            if prob[i][j-1] < 0.5 and prob[i][j] > 0.5:
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
coldJupiters = planets.loc[(planets["pl_type"] == "WJ") | (planets["pl_type"] == "CJ")]
hotJupiters = planets.loc[planets["pl_type"] == 'HJ']
seCompanions = planets.loc[planets["companion_type"] %2 == 0]
ssCompanions = planets.loc[(planets["companion_type"] % 3 == 0) | (planets["companion_type"] % 5 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]
cjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]

planetsList = [coldJupiters, hotJupiters, seCompanions, ssCompanions, hjCompanions, cjCompanions]
labels = ["Cold Jupiters", "Hot Jupiters", "SE Compnaions", "SS Companions", "HJ Companions", "CJ Companions"]

fig, ax =  plt.subplots(1,1)

for i in range(len(planetsList)):
    smax = np.logspace(np.log10(0.3),np.log10(30),49)
    cmpltMap = averageMap(np.unique(planetsList[i]["hostname"]))
    contour = findCountour(cmpltMap)
    ax.plot(smax, contour, label = labels[i])

ax.legend(frameon = False)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.3,30)
ax.set_ylim(0.3,30)
ax.set_xlabel("Semi-Major Axis (AU)")
ax.set_ylabel("Mass ($M_J$)")
plt.show()

        

    
