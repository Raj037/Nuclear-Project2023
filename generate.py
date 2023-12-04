import numpy as np
from tabulate import tabulate
import time
from HartreeFock import *
from CGgenerator import *

    
##################################################################
###                  Generate Basis States                     ###
##################################################################

# Array of parameters for gradient descent
x = []

ite = 0 # counting dummy variable

# creating a file for any given 'n' value
f = open("sdbasisGeneratedRaw.dat", 'w')
data = []
eps = []
statesDict = []

# States that should be generated
# n = [0,1,0]
# l = [2,0,2]
# j = [5,1,3]

# Start with lowest lying states
n = [0,1]
l = [2,0]
j = [5,1]

nlist = []
llist = []
jlist = []
mjlist = []
            
# Optimized Energies from literature 
# realEps = -np.asarray([4.15,4.15,8.02,8.02,3.94,3.94,7.63,7.63,3.77,3.77,6.84,6.84,2.75,2.75,3.61,3.61,8.62,8.62,11.42,11.42,7,7,9.9,9.9])
# energies for neutrons in the 0d5/2 and 1s1/2 states
realEps = -np.asarray([4.15,4.15,4.15,4.15,4.15,4.15,2.75,2.75])

# next we need to go over all possible combinations of quantum numbers, starting with 'n'
for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            # val = np.random.uniform(-2,2)
            val = -realEps[ite]
            data.append([ite,n[i],l[i],j[i],int(mj),1,val])
            eps.append(val)
            ite +=1
            # data.append([ite,n[i],l[i],j[i],int(mj),-1,val])
            # ite +=1
            statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(1)))
            # statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(-1)))
            nlist.append(n[i])
            llist.append(l[i])
            jlist.append(j[i])
            mjlist.append(mj)
            mj = mj-2
            # eps.append(val)
            # x.append(val)

# the tabulate function is slow but that should not matter for these calculations
f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

# Binding energy of 18O
converge = -np.asarray([7.767,7.767,7.767,7.767,7.767,7.767,7.767,7.767])
    
##################################################################
###           Generate Two Body Matrix Elements                ###
##################################################################


start = time.time() # I added this to help optimize the code
 
f = open("tbmeGeneratedRaw.dat", 'w')

V = np.random.rand(ite,ite,ite,ite)*2-1

for a in range(ite):
    for c in range(ite):
        for d in range(ite):
            # Antisymmetry forces these "diagonal" terms to be 0
            V[a][a][c][d] = 0
            V[a][a][d][c] = 0
            V[c][d][a][a] = 0
            V[d][c][a][a] = 0
            
for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                    
                        # Randomizing starting guess by pulling from a normal distribution
                        # centered around formula in class with a width of 12% (2 sigma?)
                        if a < j[0]+1 and b < j[0]+1 and c < j[0]+1 and d < j[0]+1:
                            ME = np.random.normal(3.91,3.91*0.12,1)
                        elif a > j[0] and b > j[0] and c > j[0] and d > j[0]:
                            ME = np.random.normal(3.65,3.65*0.12,1)
                        else:
                            ME = np.random.normal(2.68,2.68*0.12,1)
                            
                        # Normalizing using Clebsch-Gordon Coefficients (Take a look at CGgenerator.py)
                        CG = 0
                        Mp = abs(int(mjlist[c]/2+mjlist[d]/2))
                        Jp = int(jlist[c]/2+jlist[d]/2)
                        for JM in range(Mp,Jp+1):
                            CG += CGcoeff(jlist[c]/2,jlist[d]/2,mjlist[c]/2,mjlist[d]/2,JM,Mp)
                            
                        ME = CG*ME
                        
                        # Add to parameter array
                        x.append(ME[0])
                        
                        # Antisymmetry
                        V[a][b][c][d] = ME[0]
                        V[b][a][d][c] = ME[0]
                        V[b][a][c][d] = -ME[0]
                        V[a][b][d][c] = -ME[0]
                        
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(ME[0]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(-ME[0]) + '\n')
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(-ME[0]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(ME[0]) + '\n')
 						
f.close()
end = time.time()
 
#print(end-start)

    
##################################################################
###        Optimize Energies and Two Body Matrix Elements      ###
##################################################################

step = 1e-3
# rate = 1e-3
rate = np.zeros(len(x))
A = 2
Error = 1e20
dE = 1e-3
stepN = 0
maxSteps = 400
max_iter = 20

# The learning rate for greatest descent
# This gets the magnitude of N, ie returns x in 10^x
def getRate(N):
    return np.floor(np.log10(abs(N)))

# Converts array of parameters into an array of hamiltonian energies and tbm
def x2epsV(vec_x):
    energy = np.zeros(ite)
    tbme = np.zeros([ite,ite,ite,ite])
    
    ind = 0
    
# =============================================================================
#     for i in range(int(ite/2)):
#         energy[2*i] = vec_x[ind]
#         energy[2*i+1] = vec_x[ind]
#         
#         ind += 1
# =============================================================================
        
    for a in range(ite):
        for b in range(ite):
            if (a < b):
                for c in range(ite):
                    for d in range(ite):
                        if (c < d):
                            
                            ME = vec_x[ind]
                            ind += 1
                                
                            tbme[a][b][c][d] = ME
                            tbme[b][a][d][c] = ME
                            tbme[b][a][c][d] = -ME
                            tbme[a][b][d][c] = -ME
                            
    return energy, tbme


# The functions to minimize is the error function, fitting least squares to the 
# binding energy of 18O
def calcError(E):
    return np.sum(np.square(E[:A]-converge[:A]))

oldError = 0
diffError = 1
        
# Find the gradient for each element

while Error > dE and stepN < maxSteps:

    df = np.zeros(len(x))

    for i in range(len(x)):
        
        xPlus = x.copy()
        xMinus = x.copy()
        
        # For each parameter, go both directions
        xPlus[i] += step
        xMinus[i] -= step
        
        # Parameters to matrices
        skip, VPlus = x2epsV(xPlus)
        skip, VMinus = x2epsV(xMinus)
        
        # Hartree Fock to get the energies
        epsNewPlus, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VPlus)    
        epsNewMinus, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VMinus)
        
        # Calculate the error
        fplus = calcError(epsNewPlus[-1])
        fminus = calcError(epsNewMinus[-1])
        
        # Get the slope/gradient
        df[i] = (fplus-fminus)/(2*step)
        
        # Dynamic learning rate, here we change parameters by 10^-3
        rate[i] = 10**(-(getRate(df[i])+3))
        
        # To speed up the process, at the beginning, change parameters more,
        # and near the end, get better precision
        if stepN < 10:
            rate[i] = rate[i]*3
        if stepN > 15:
            rate[i] = rate[i]/2
        
    # Change the parameters by going in the negative gradient direction
    x = x-rate*df
    
    # Get new tbme
    skip, VN = x2epsV(x)
    
    # Used to calculate the current error
    epsNew, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VN)
    
    oldError = Error
    Error = calcError(epsNew[-1])
    
    diffError = abs(Error-oldError)
    
    print(Error)
    print(epsNew[-1][:A])
    
    stepN += 1
    
# Prepare and print the final elements
skip, finalV = x2epsV(x)

finalEps = epsNew.copy()

f = open("sdbasisGeneratedOptimized.dat","w")
dataOpt = []
k = 0

for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),1,finalEps[k]])
            k += 1
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),-1,finalEps[k]])
            k += 1
            mj = mj-2
            

f.write(tabulate(dataOpt, tablefmt="plain",showindex=False))
f.close()
    
f = open("tbmeGeneratedOptimized.dat","w")

for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(finalV[a][b][c][d]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(finalV[b][a][c][d]) + '\n')
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(finalV[a][b][d][c]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(finalV[b][a][d][c]) + '\n')
                        
f.close()