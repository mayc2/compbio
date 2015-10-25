import numpy
import scipy
import math
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO

#parse input file into record data structure
record = SeqIO.read(open("sequence.fasta"),"fasta")

#window sizes
w = [500, 1000, int(len(record.seq)/20)]

#find number of points to plot for each window size (x values)
length = len(record.seq)
indices = []
for window in w:
	indices.append(int(math.ceil(length/window)))

#initializing arrays that hold freq at each window size
a_arr = []
c_arr = []
t_arr = []
g_arr = []
at_arr = []
cg_arr = []
x_values = []

#for each window, add a list
for i in range(0,len(w)):

	a_arr.append([])
	c_arr.append([])
	t_arr.append([])
	g_arr.append([])
	at_arr.append([])
	cg_arr.append([])
	x_values.append([])

	#index of current position in sliding window
	start = 0

	#for each window in the seq, find the freq 
	for x in range(0, indices[i]):
		x_values[i].append(start)
		end = start + (w[i]-1)
		if end > (length - 1):
			end = (length - 1)
		
		a_arr[i].append(float(record.seq[start:end].count("A"))/float(w[i]))
		c_arr[i].append(float(record.seq[start:end].count("C"))/float(w[i]))
		t_arr[i].append(float(record.seq[start:end].count("T"))/float(w[i]))
		g_arr[i].append(float(record.seq[start:end].count("G"))/float(w[i]))
		
		at_arr[i].append((float(record.seq[start:end].count("A"))+float(record.seq[start:end].count("T")))/float(w[i]))
		cg_arr[i].append((float(record.seq[start:end].count("C"))+float(record.seq[start:end].count("G")))/float(w[i]))
		start = end + 1

titles = []
titles2 = []
for i in range(0,len(w)):
	titles.append("Frequency of A,C,T,G Nucleotides at Window Size " + str(w[i]))
	titles2.append("Frequency of A-T, C-G Nucleotides at Window Size " + str(w[i]))

#plot windows
for i in range(0,len(w)):
	plt.plot(x_values[i],a_arr[i],label="A Nucleotide")
	plt.plot(x_values[i],c_arr[i],label="C Nucleotide")
	plt.plot(x_values[i],t_arr[i],label="T Nucleotide")
	plt.plot(x_values[i],g_arr[i],label="G Nucleotide")
	plt.title(titles[i])
	plt.ylabel("Frequency")
	plt.xlabel("Index")
	plt.ylim(0.0,1.0)
	plt.legend()
	plt.show()

for i in range(0,len(w)):
	plt.plot(x_values[i],at_arr[i],label="A-T Nucleotides")
	plt.plot(x_values[i],cg_arr[i],label="C-G Nucleotides")
	plt.title(titles2[i])
	plt.ylabel("Frequency")
	plt.xlabel("Index")
	plt.ylim(0.0,1.0)
	plt.legend()
	plt.show()