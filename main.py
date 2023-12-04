import numpy as np
from HartreeFock import *

##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################

# Data is stored in dictionaries for explicit labeling

statefile = open('spbasis.dat','r')

statesDict = []

for x in statefile:
    temp1,temp2,temp3,temp4,temp5,temp6 = x.split()
    statesDict.append(dict(n=int(temp2),l=int(temp3),j=int(temp4),mj=int(temp5),t3=int(temp6)))
    
lenStates = len(statesDict)
    
Vfile = open('tbme.dat','r')

VDict = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    VDict.append(dict(a=int(temp1)-1,b=int(temp2)-1,c=int(temp3)-1,d=int(temp4)-1,val=float(temp5)))
    
# Calculate N=2n+l for all states, not needed for sd-basis
for s in statesDict:
    s.update({"N":2*s["n"]+s["l"]})    
    
##################################################################
###                 INITIALIZING MATRICES                      ###
##################################################################
    
# Setting up matrices
    
V = np.zeros((lenStates,lenStates,lenStates,lenStates))

for x in VDict:
    V[x["a"]][x["b"]][x["c"]][x["d"]] = x["val"]

# Energies of Harmonic Oscillator
eps = np.zeros(lenStates)

hw = 10 # HO quanta, chosen in MeV
for i in range(lenStates):
    eps[i] = hw*(statesDict[i]["N"]+3/2)
    
##################################################################
###                    Run Hartree Fock                        ###
##################################################################
    
# Running Hartree Fock
# Function is stored in HartreeFock.py
# Inputs are HartreeFock( A , Smallest Difference in Energy Wanted ,
#                           Max Number of iterations requests , 
#                           Basis Stored as a Dict , Energies in an Array,
#                           and Two Body Density Matrix as 4d Array)
# Outputs are bindingEnergies, a 2d array where the first index is the iteration number
#                           and the second index is the state index
# To access the final iteration, use bindingEnergies[-1] as the array
# Expected negEnergies = A

# A is number of nucleons
A = 4

bindingEnergies, negEnergies = HartreeFock(A,1e-5,1000,statesDict,eps,V)

# Writing out final energies
# orbit defined by l quantum number
orbit = ["s","p","d","f","g","h"]
f = open("BindingEnergiesHO.txt","w")

for i in range(lenStates):
    # prints out in orbital notation nl^pi,nu (j)
    if statesDict[i]["t3"] == 1: nuc = "n"
    else: nuc = "p"
    f.write(str(statesDict[i]["n"])+orbit[statesDict[i]["l"]]+"^"+nuc+"("+str(statesDict[i]["j"])+"/2) " + str(round(bindingEnergies[-1][i],4))+ " MeV\n")
    
B = np.zeros(len(bindingEnergies))

# Calculates the binding energy per nucleon and stores it in an array
for i in range(len(bindingEnergies)):
    B[i] = sum(bindingEnergies[i][:negEnergies])/negEnergies
        

f.close()
