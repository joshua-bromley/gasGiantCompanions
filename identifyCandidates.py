import pandas as pd
import numpy as np

seenCandidates = pd.read_csv("./data/sampleCandidates.csv", skiprows = 100)
candidates = pd.read_csv("./data/sampleCandidatesSP.csv", skiprows=99)
sample = pd.read_csv("./data/gasGiantDataComplete2.csv")

seenSystems = np.unique(seenCandidates["hostname"])
candidateSystems = np.unique(candidates["hostname"])
sampleSystems = np.unique(sample["hostname"])

remainingCandidates = []
inSample = []

for i in range(len(candidateSystems)):
    if not candidateSystems[i] in sampleSystems and not candidateSystems[i] in seenSystems:
        remainingCandidates.append(candidateSystems[i])
    else:
        inSample.append(candidateSystems[i])

np.savetxt("./sampleCandidatesSP.txt", remainingCandidates, fmt = "%s")