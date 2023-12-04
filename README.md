# NUCPROJ

## Directory Structure
* main.py
  * Reads in spbasis.dat and tbme.dat, then uses Hartree-Fock to calculate the single particle energies for a given A.
* HartreeFock.py
  * The main function for Hartree-Fock, can be imported to any file for general use, ex. main.py.
* generate.py
  * Generates  $0d_{5/2}$  and  $1s_{1/2}$  states, randomly generates two body matrix elements using a normal distribution, and then optimizes the tbme using greatest descent.
* CGgenerator.py
  * Generates the Clebsch-Gordon coefficients to normalize the tbme.
 
## HartreeFock.py
Given A (number of nucleons), epsilon (minimum value for convergence), max_iter (maximum number of iterations), a dictionary of states, h_0 energies, and two-body matrix elements, performs the Hartree-Fock method. 

Starts by generating the $C_{ij}$ matrices, which can be random (commented out). We instead, for consistency during optimization, start with an identity matrix. Starting the loop, uses $C_{ij}$ to calculate the $\rho_{\delta\gamma}$ matrix. Then, using the tbme and $\rho$, calculate the HF potential, which is combined with the $h_0$ energies to get the total HF Hamiltonian. Using the HF Hamiltonian, uses linear algebra to get the eigenvalues and eigenvectors. Use the eigenvectors as the new $C$ in the next loop, and use the eigenvalues to calculate the change in energy.

## main.py
Starts by importing spbasis.dat and tbme.dat as dictionaries to make reading easier, then calculates the HO energies from n and l. Passes everything through HartreeFock and returns a list of all iterations. At the end, we calculate B/A for all iterations.

## CGgenerator.py
Checks for negative factorials, and there are none, calculates the CG coefficients.

## generate.py
Generates all possible  $0d_{5/2}$  and   $1s_{1/2}$  states, ie. finds all possible $m_j$ and $t_3$ for level. We ignore  $0d_{3/2}$  since we start with the lowest shell. We also ignore all protons (we only set $t_3$ = 1) for simplicity. The one body matrix elements can be set using physical values. For example: 
$$e_{0d_{5/2}}=E(^{17}O,\frac{5}{2}^+)-E(^{16}O,O^+)$$
The tbme are randomly generated from a normal distribution using a mean of some similar physical values and a width of 12%. For example, the mean for the $0d_{5/2}$ can be found as:
$$<(0d5/2)^2|\hat{V}|(0d5/2)^2)>=E(^{18}O,O^+)-2e_{0d5/2}-E(^{16}O,O^+)$$
From here, we start optimizing these tbme while fixing the one body matrix elements. The function we try to minimize is (least squares minimization):
$$f(\vec{x}) = \sum_{i}(x_i - x_{data})^2$$
where $x_data$ is the experimental single particle energies and $x$ is the calculated single particle energies, and $i$ ranges from $1$ to $A$. 

The optimization procedure starts with finding the gradient of our function. This is done one at a time for each parameter in our array. Then we change the parameters by going against the gradient using a learning rate, which we dynamically set to speed up the process. We set learning rate such that we change the parameters by a magnitude of $10^{-3}$. The loop ends once we find an error below some threshold.
