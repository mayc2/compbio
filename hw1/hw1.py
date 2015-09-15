from Bio import SeqIO
import numpy
import scipy
import matplotlib.pyplot as plt

for n in SeqIO.parse("sequence.fasta","fasta"):
	next

#test
#print(n.seq)
arr = []
arr.append(n.seq.count("A"))
arr.append(n.seq.count("C"))
arr.append(n.seq.count("T"))
arr.append(n.seq.count("G"))

at = arr[0] + arr[2]
cg = arr[1] + arr[3]

#print counts of each
#print(arr, at, cg)

#window sizes
w = [500, 1000, (len(n.seq)/20)]

#print window sizes
#print(w)

#find number of points to plot for each window size (x values)
length = len(n.seq)
print(length)
indices = []
for window in w:
	indices.append(length/window)
print(indices)

#plot windows
#plt.plot(12)
#plt.ylabel("Count")
#plt.show()