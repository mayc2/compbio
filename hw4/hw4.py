import sys
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import random
from Bio import SeqIO

def parse_motifs(filename):
	records = list(SeqIO.parse(open(filename),"fasta"))
	motifs = []
	start_indexes = []
	sequences = []

	for sequence in records:
		temp = sequence.description.split()
		motifs.append(temp[2])
		start_indexes.append(int(temp[1]))
		sequences.append(sequence.seq)
	return sequences,motifs,start_indexes

def motif_model(motifs,exclude):
	#pseudo count = 1.0
	answer = g_init(len(motifs[0]))
	for j in range(len(motifs[0])):
		A_count = 0.0
		C_count = 0.0
		T_count = 0.0
		G_count = 0.0
		for i in range(len(motifs)):
			if (i != exclude and j <= 92):
				if(motifs[i][j] == 'A'):
					A_count += 1.0
				if(motifs[i][j] == 'C'):
					C_count += 1.0
				if(motifs[i][j] == 'T'):
					T_count += 1.0
				if(motifs[i][j] == 'G'):
					G_count += 1.0
		answer['A'][j] = (A_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
		answer['C'][j] = (C_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
		answer['T'][j] = (T_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
		answer['G'][j] = (G_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
	return answer

def background_model(start_indexes,sequences,exclude):
	answer = {}
	for x in ['A','C','T','G']:
		answer[x] = 0.0

	A_count = 0.0
	C_count = 0.0
	T_count = 0.0
	G_count = 0.0

	#for each position in the sequence
	for j in range(len(sequences[0])):
		#for each sequence
		for i in range(len(sequences)):
			ignore = int(start_indexes[i])-1
			if (i != exclude):
				if ((j < ignore or j >= (ignore + 8)) and (j <= (len(sequences[0])-8))):
						if(sequences[i][j] == 'A'):
							A_count += 1.0
						if(sequences[i][j] == 'C'):
							C_count += 1.0
						if(sequences[i][j] == 'T'):
							T_count += 1.0
						if(sequences[i][j] == 'G'):
							G_count += 1.0
	answer['A'] = (A_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
	answer['C'] = (C_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
	answer['T'] = (T_count + 1.0) / (A_count+C_count+T_count+G_count+4.0) 
	answer['G'] = (G_count + 1.0) / (A_count+C_count+T_count+G_count+4.0)
	return answer

def find_pm(sequence,Theta_m):
	answer = []
	for j in range(len(sequence)-7):
		temp = 1.0
		for i in range(j,j+8):
			temp *= Theta_m[sequence[i]][i-j]
		answer.append(temp) 
	return answer

def find_pb(sequence,Theta_B):
	answer = []
	for j in range(len(sequence)-7):
		answer.append(Theta_B[sequence[j]]) 
	return answer

def cal_ratios(Probability_m,Probability_B):
	answer = []
	for i in range(len(Probability_m)):
		answer.append(Probability_m[i]/Probability_B[i])
	return answer

def normalize(ratios):
	total = sum(ratios[i] for i in range(len(ratios)))
	for i in range(len(ratios)):
		ratios[i] /= total
	return ratios

def gibbs(sequences,motifs,start_indexes,burnin):

	w = len(motifs[0])
	L = len(sequences[0])

	#initialize estimated motif start sites
	A = []
	for i in range(len(motifs)):
		A.append(random.randrange(0,L-w))

	#burnin iteration
	for _ in range(burnin):
		for i in range(len(motifs)):
			a_temp = 0.0

			#Motif Model Calculation
			Theta_m = motif_model(motifs,i)

			#Background Model Calculation
			Theta_B = background_model(A,sequences,i)

			#Probability from Motif Model
			Probability_m = find_pm(sequences[i],Theta_m)

			#Probability from Background Model
			Probability_B = find_pb(sequences[i],Theta_B)

			#ratio calculation
			ratios = cal_ratios(Probability_m,Probability_B)

			#estimate new motif start location
			ratios = normalize(ratios)
			a_temp = np.random.choice(93,p=ratios)
			A[i] = a_temp
	Burn = A

	# sampling iteration
	C = c_init()
	for _ in range(burnin*2):
		for i in range(len(motifs)):
			a_temp = 0.0

			#Motif Model Calculation
			Theta_m = motif_model(motifs,i)

			#Background Model Calculation
			Theta_B = background_model(A,sequences,i)

			#Probability from Motif Model
			Probability_m = find_pm(sequences[i],Theta_m)

			#Probability from Background Model
			Probability_B = find_pb(sequences[i],Theta_B)

			#ratio calculation
			ratios = cal_ratios(Probability_m,Probability_B)

			#estimate new motif start location
			ratios = normalize(ratios)
			a_temp = np.random.choice(93,p=ratios)
			A[i] = a_temp

			#update Count
			C[i][a_temp] += 1.0
	return Burn,A,C,Theta_m,Theta_B

def c_init():
	answer = []
	for i in range(10):
		temp = []
		for j in range(100):
			temp.append(0.0)
		answer.append(temp)
	return answer

def g_init(w):
	answer = {}
	for i in ['A','C','T','G']:
		answer[i] = []
		for j in range(w):
			answer[i].append(0.0)
	return answer

def normalize_C(C,count):
	answer = []
	for i in range(len(C)):
		temp = []
		for j in range(len(C[i])):
			temp.append(float(C[i][j])/float(count))
		answer.append(temp)
	return answer

def main():

	#parse the sequences,motifs and start indexes from a record into lists
	sequences,motifs,start_indexes = parse_motifs("test.fasta")
	
	#set burnin and iterate through the gibbs sampler 5 times
	burnin = 1000
	Burn = []
	g_indexes = []
	C = []
	m = []
	B = []
	for _ in range(5):
		t1,t2,t3,t4,t5 = gibbs(sequences,motifs,start_indexes,burnin)
		Burn.append(t1)
		g_indexes.append(t2)
		C.append(t3)
		m.append(t4)
		B.append(t5)

	#sum counts and normalize
	Counts = c_init()
	for iteration in range(len(C)):
		for i in range(len(C[iteration])):
			for j in range(len(C[iteration][i])):
				Counts[i][j] += C[iteration][i][j]
	Counts = normalize_C(Counts,5*2*burnin)

	#average motif model
	Theta_m = g_init(10)
	for i in ['A','C','T','G']:
		for j in range(len(m[0][i])):
			temp = 0.0
			for iteration in range(len(m)):
				temp += m[iteration][i][j]
			Theta_m[i][j] = temp/5.0
	print("\nAverage Motif Model")
	for i in ['A','C','T','G']:
		temp = str(Theta_m[i][0])
		for j in range(1,8):
			temp = temp + " " + str(Theta_m[i][j])
		print("\n" + i + ": " + temp)

	#average of background model
	Theta_B = {}
	for i in ['A','C','T','G']:
		Theta_B[i] = 0.0
	for i in ['A','C','T','G']:
		temp = 0.0
		for iteration in range(len(B)):
			temp += B[iteration][i]
		Theta_B[i] = temp/5.0
	print("\nAverage Background Model")
	for i in ['A','C','T','G']:
		print(i + " " + str(Theta_B[i]))

	#print Motif Locations, Sequences and Frequencies for each sequence
	print("\nMotif Locations, Sequences and Frequencies")
	for seq in range(len(sequences)):
		print("Seq " + str(seq+1) + ": " + str(start_indexes[seq]) + " " + motifs[seq] + " " + str(Counts[seq][start_indexes[seq]-1]) )

	#display graphs for the probabilities of each sequence
	for i in range(len(Counts)):
		plt.bar(np.arange(1,101),Counts[i],label="Count for Seq "+str(i+1))
		plt.title("Sampling Probability for Sequence "+str(i+1))
		plt.ylabel("Frequency")
		plt.xlabel("Position")
		plt.xlim(0,100)
		plt.ylim(0,1)
		plt.show()

	return 0

if __name__ == "__main__":
	sys.exit(main())