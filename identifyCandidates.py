import pandas as pd
import numpy as np

candidates = pd.read_csv("./data/sampleCandidates.csv", skiprows=100)
sample = pd.read_csv("./data/gasGiantDataComplete2.csv")

candidateSystems = np.unique(candidates["hostname"])
sampleSystems = np.unique(sample["hostname"])

remainingCandidates = []

for i in range(len(candidateSystems)):
    if not candidateSystems[i] in sampleSystems:
        remainingCandidates.append(candidateSystems[i])

np.savetxt("./sampleCandidates.txt", remainingCandidates, fmt = "%s")