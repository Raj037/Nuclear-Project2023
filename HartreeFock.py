import numpy as np

def HartreeFock(A,sat,max_iter,states,eps,V):
    
    numStates = len(eps)
    
    dE = 1
    it_count = 0
    ener_old = np.zeros(numStates)
    energies = np.zeros(numStates)
    orbit = ["s","p","d","f","g","h"]
    
    allEnergies = []
    
    # initial C matrix is guessed at random
    # C = np.random.rand(numStates,numStates)    
    C = np.eye(numStates)
    
    # looping until the difference is minimal
    while ((dE > sat) and (it_count < max_iter)):
        # Calculate the one body density matrix
        rho = np.zeros((numStates,numStates))

        # Calculating the one body density matrix, rho
        for delta in range(numStates):
            for gamma in range(numStates):
                for j in range(A):
                    # Make sure that j does not exceed the size of C or V when choosing A
                    if j>=numStates: continue;
                    rho[delta][gamma] += C[delta][j]*C[gamma][j]
        
        # Calculate HF Potential
        U = np.zeros((numStates,numStates))
        
        for eta in range(numStates):
            for beta in range(numStates):
                for delta in range(numStates):
                    for gamma in range(numStates):
                        U[eta][beta] += rho[delta][gamma]*V[eta][gamma][beta][delta]
                        
        # For initial <eta|h_0|gamma>, use diagonal matrix with energies
        E = np.diag(eps)
                
        # Calculate Hartree-Fock Hamiltonian
        HF = E+U
        
        # Calculate new C matrix and energies, store previous energies
        ener_old = energies
        # Eigenvalue problem solved by linalg package in numpy
        energies, C = np.linalg.eigh(HF)
        
        # Calculate energy difference
        dE = np.sum(np.abs(energies-ener_old))/numStates
        it_count += 1
        
# =============================================================================
#         print("Iteration " + str(it_count))
#         for i in range(numStates):
#             # prints out in orbital notation nl^pi,nu (j)
#             if states[i]["t3"] == 1: nuc = "n"
#             else: nuc = "p"
#             print(str(states[i]["n"])+orbit[states[i]["l"]]+"^"+nuc+"("+str(states[i]["j"])+"/2) " + str(round(energies[i],4))+ " MeV")
#             
#         print("\n")
# =============================================================================
        
        allEnergies.append(energies)
        
    countNeg = 0
    for i in range(len(allEnergies[-1])):
        if allEnergies[-1][i] < 0:
            countNeg += 1
        
    return allEnergies, countNeg
