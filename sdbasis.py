import numpy as np
from tabulate import tabulate

# First we need 2to decide the size of the model space by selecting an 'n' value
n_max = int(input('Give n(max like 3 tbh): '))

# we are working in sd so l = 0 or 2
l = [0,2]

ite = 1 # counting dummy variable

# creating a file for any given 'n' value
f = open("sdbasis_nmax" + str(n_max) + ".txt", 'w')
data = []

# next we need to go over all possible combinations of quantum numbers, starting with 'n'
for i in range(len(l)): # there is only l=0 and l=2 in this case but could be adjusted for any number
	n = 0
	while n <= n_max: 
		j = l[i] + 1 / 2 # definition of total angular momentum
		mj_min = -j # mj goes from j to -j in increments of 1
		mj = j
		while mj >= mj_min:
			for c in range(2): # we use if statements to set the isospin
				if c == 0:
					data.append([ite,n,l[i],int(l[i]+ 1 % 2),int(2*mj),1])
					ite +=1
				if c == 1:
					data.append([ite,n,l[i],int(l[i] + 1 % 2),int(2* mj),-1])
					ite +=1
			mj = mj-1
		n += 1	

# the tabulate function is slow but that should not matter for these calculations
f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

