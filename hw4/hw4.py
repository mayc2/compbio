import sys
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import random
from Bio import SeqIO

def parse(filename):

	answer = list(SeqIO.parse(open(filename),"fasta"))
	return answer

def parse_motifs(records):
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
	# print (A)

	#burnin iteration
	for _ in range(burnin):
		for i in range(len(motifs)):
			a_temp = 0.0

			#Motif Model Calculation
			Theta_m = motif_model(motifs,i)
			# print(Theta_m)

			#Background Model Calculation
			Theta_B = background_model(start_indexes,sequences,i)
			# print(Theta_B)

			#Probability from Motif Model
			Probability_m = find_pm(sequences[i],Theta_m)
			# print(Probability_m)

			#Probability from Background Model
			Probability_B = find_pb(sequences[i],Theta_B)
			# print(Probability_B)

			#ratio calculation
			ratios = cal_ratios(Probability_m,Probability_B)
			# print(ratios)

			#estimate new motif start location
			ratios = normalize(ratios)
			a_temp = np.random.choice(93,p=ratios)
			A[i] = a_temp
			# print(a_temp)
	Burn = A
	# sampling iteration
	C = c_init()
	for _ in range(burnin*2):
		for i in range(len(motifs)):
			a_temp = 0.0

			#Motif Model Calculation
			Theta_m = motif_model(motifs,i)
			# print(Theta_m)

			#Background Model Calculation
			Theta_B = background_model(start_indexes,sequences,i)
			# print(Theta_B)

			#Probability from Motif Model
			Probability_m = find_pm(sequences[i],Theta_m)
			# print(Probability_m)

			#Probability from Background Model
			Probability_B = find_pb(sequences[i],Theta_B)
			# print(Probability_B)

			#ratio calculation
			ratios = cal_ratios(Probability_m,Probability_B)
			# print(ratios)

			#estimate new motif start location
			ratios = normalize(ratios)
			a_temp = np.random.choice(93,p=ratios)
			A[i] = a_temp
			C[i][a_temp] += 1.0
			# print(a_temp)
	return Burn,A,C

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

"""
When to quit?
posterior alignment probability, 
quit when it plateaus, 
report MAP(max posteriori probability)
give point estimate of models and site positions

i = sequence
j = position

Motif Model
Theta_m = (n[i][j] + a[i][j]) / sum(n[i][l] + a[i][l] for all l)

Background Model
Theta_B = (n[B][j] + a[B][j]) / sum(n[B][l] + a[B][j] for all l) 

Sample position based on probability ratio of each position
P(X[i][j]...X[i][j+w-1] | Theta_m) / P(X[i][j]...X[i][j+w-1] | Theta_B)

"""
def normalize_C(C,count):
	answer = []
	for i in range(len(C)):
		temp = []
		for j in range(len(C[i])):
			temp.append(float(C[i][j])/float(count))
		answer.append(temp)
	return answer

def adjust_index(Indexes):
	answer = []
	for index in Indexes:
		answer.append(index + 1)
	return answer

def main():

	records = parse("test.fasta")
	# print (records)

	sequences,motifs,start_indexes = parse_motifs(records)
	# print (motifs)
	
	burnin = 100
	print("Calculating for a burn-in:", burnin)
	Burn,g_indexes,C = gibbs(sequences,motifs,start_indexes,burnin)
	Burn = adjust_index(Burn)
	g_indexes = adjust_index(g_indexes)
	print("Actual Start Indexes:")
	print(start_indexes)
	print("Estimated Burn Start Indexes:")
	print(Burn)
	print("Estimated Sample Start Indexes:")
	print(g_indexes)

	C = normalize_C(C, burnin*2)

	for i in range(len(C)):
		plt.bar(np.arange(1,101),C[i],label="Count for Seq "+str(i))
		plt.title("Sampling Probability for Sequence "+str(i+1))
		plt.ylabel("Frequency")
		plt.xlabel("Position")
		plt.xlim(0,100)
		plt.show()

	# burnin = 1000
	# print("\nCalculating for a burn-in:", burnin)
	# Burn,g_indexes,C1 = gibbs(sequences,motifs,start_indexes,burnin)
	# Burn = adjust_index(Burn)
	# g_indexes = adjust_index(g_indexes)
	# print("Actual Start Indexes:")
	# print(start_indexes)
	# print("Estimated Burn Start Indexes:")
	# print(Burn)
	# print("Estimated Sample Start Indexes:")
	# print(g_indexes)

	# C1 = normalize_C(C1, burnin*2)

	# for i in range(len(C1)):
	# 	plt.bar(np.arange(1,101),C1[i],label="Count for Seq "+str(i))
	# 	plt.title("Sampling Probability for Sequence "+str(i+1))
	# 	plt.ylabel("Frequency")
	# 	plt.xlabel("Position")
	# 	plt.xlim(0,100)
	# 	plt.show()

	return 0

if __name__ == "__main__":
	sys.exit(main())