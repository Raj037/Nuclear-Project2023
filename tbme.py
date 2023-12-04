import numpy as np
import time

start = time.time() # I added this to help optimize the code

# we need to start by determining the amount of quanta each value (a,b,c,d) has
# this is determined by the model space size
# from the sdbasis, q= # of states
q = int(input("Enter quanta number: "))
f = open("tbme_" + str(q)+ ".txt", 'w')
table = []

# quadruple for loops go brr
for a in range(1,q+1):
	for b in range(1,q+1):
		if (a != b and a < b): # first if statement that determines which values of a and b are important - helps run faster
			for c in range(1,q+1):
				for d in range(1, q+1):
					if (c !=d and c < d): # second if statement determines which c and d are important
						val = np.random.uniform(-2,2) # ussing antisymmetry, we can assign four states to the same random value
						f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(val) + '\n')
						f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(-val) + '\n')
						f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(-val) + '\n')
						f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(val) + '\n')
						

f.close()
end = time.time()

print(end-start)
