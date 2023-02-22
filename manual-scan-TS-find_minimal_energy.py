# Author: disa
# Description:This program find every energy of each cycle that gaussian made analyzing a molecule, then sort the energies and extract the lower point of energy (lower value) for all the files and generate a the energies.dat file. Next it takes the txt file to sor the energy from the low to high per file generating the minimal_energies.txt 
# Instructions: python3 find_minimal_energy-py

import numpy as np
import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

argv = sys.argv
path_input = argv[1]
path_output = f"out_{path_input}"
os.system(f"mkdir {path_output}")
os.system(f"mv {path_input} {path_output}")

# List log files.
logFilesList = []
for fileName in os.listdir(path_output+path_input):
    if fileName.endswith('.log'): logFilesList.append(fileName)

# Find and list energies.
fileNumbers, minEnergies = [], []
for logFileName in logFilesList:
    with open(path_output+path_input+logFileName,'r') as logFile: fileLines = logFile.readlines()
    energies = []
    for line in fileLines:
        findEnergy = re.search(r'E\(RB3LYP\)\s+=\s+(-?\d+\.?\d*)\s+',line)
        if findEnergy: 
            energies.append(float(findEnergy.group(1)))
    if len(energies) == 0: 
        print(f'Warning: E(RB3LYP) not found in file {logFileName}')
        continue
    energies.sort() 
    minEnergies.append(energies[1])
    fileNumbers.append(int(re.sub(r'.+?(\d+)\.log',r'\1',logFileName)))
fileNumbers = np.array(fileNumbers)
minEnergies = np.array(minEnergies) #Hartree

# Create data frame.
header = ['File number','Energy[Hartree]']
data = np.vstack((fileNumbers,minEnergies)).T
dataFrame = pd.DataFrame(data, columns=header)
dataFrame['File number'] = dataFrame['File number'].map(int)
dataFrame.sort_values(['Energy[Hartree]'],inplace=True,ignore_index=True)
lowest_energy = dataFrame['Energy[Hartree]'].min()
dataFrame['DeltaEnergy[Hartree]'] = dataFrame['Energy[Hartree]'] - lowest_energy
dataFrame['Energy[kcal/mol]'] = dataFrame['Energy[Hartree]']*627.50947406
dataFrame['DeltaEnergy[kcal/mol]'] = dataFrame['DeltaEnergy[Hartree]']*627.50947406

# Print in screen average energy and lowest energy log file
lowest_energy_row = dataFrame[['Energy[Hartree]','Energy[kcal/mol]']][dataFrame['Energy[Hartree]'] == lowest_energy]
average_energy_Hartrees = dataFrame['Energy[Hartree]'].mean()
average_energy_kcalmol = dataFrame['Energy[kcal/mol]'].mean()
print(f'Minimal energy \n'
      f'{lowest_energy_row}, \n'
      f'Average energy {average_energy_Hartrees} Hartrees, \n'
      f'Average energy {average_energy_kcalmol} kcal/mol')

#Plot energies for MANUAL SCAN TS
dataFrame.sort_values(['File number'],inplace=True,ignore_index=True)
x = dataFrame['File number'].values
y = dataFrame['Energy[Hartree]'].values
mymodel = np.poly1d(np.polyfit(x, y, 10))
fig,ax = plt.subplots(1)
x = np.linspace(np.amin(x),np.amax(x),100)
ax.plot(x,mymodel(x))
dataFrame.plot(x = 'File number', y = 'Energy[Hartree]', kind = 'scatter', ax=ax, xlabel = 'distance[100 A]')
plt.tight_layout()
#plt.show()
plt.savefig(f'{path_output}scan-TS.pdf')

# Save energies.
dataFrame.to_csv(f'{path_output}minimal_energies_DFT.dat',sep='\t',index=False)
