import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import optimize as opt

planets = pd.read_csv("./data/gasGiantData.csv")

coldJupiters = planets.loc[planets["pl_type"] == "CJ"]
wcjCompanions = planets.loc[(planets["companion_type"] % 11 == 0) | (planets["companion_type"] % 13 == 0)]
hjCompanions = planets.loc[planets["companion_type"] % 7 == 0]

gasGiantsList = [coldJupiters, wcjCompanions, hjCompanions]

bins = np.logspace(np.log10(0.5),np.log10(25),11)
binCentres = np.sqrt(np.multiply(bins[:10],bins[1:]))
logBinCentres = np.log10(binCentres)
binCounts = np.zeros((len(gasGiantsList),10))

for k in range(len(gasGiantsList)):
    for i in range(len(gasGiantsList[k])):
        smax = gasGiantsList[k].iloc[i]["pl_bmassj"]
        for j in range(len(binCounts[k])):
            if smax > bins[j] and smax <= bins[j+1]:
                binCounts[k][j] += 1
                break
    
binErrors = np.sqrt(binCounts)

for i in range(len(gasGiantsList)):
    binCounts[i] /= len(gasGiantsList[i])
    binErrors[i] /= len(gasGiantsList[i])

def linearModel(x, m, b):
    return m*x + b


params = np.zeros((len(gasGiantsList),2))
for i in range(len(gasGiantsList)):
    param, cov = opt.curve_fit(linearModel, logBinCentres, binCounts[i], sigma = binErrors[i])
    print(param)
    params[i] = param

x = np.log10(np.logspace(0,np.log10(20),100))
fig, ax = plt.subplots(1,1)
for i in range(len(gasGiantsList)):
    ax.errorbar(logBinCentres, binCounts[i], yerr = binErrors[i], ls = "", marker = "o")
    ax.plot(x, linearModel(x, *params[i]))


plt.show()




