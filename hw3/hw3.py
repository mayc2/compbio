import sys
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import random
from Bio import SeqIO

def parse(inFile):
	answer = list(SeqIO.parse(open(inFile), "fasta"))
	# print (answer)
	return answer

#use random to findi initial values
def emission():
	answer = []
	for _ in range(2):
		temp = random.random()
		temp2 = random.random()
		temp3 = random.random()
		temp4 = random.random()
		n = temp+temp2+temp3+temp4
		temp /= n
		temp2 /= n
		temp3 /= n
		temp4 /= n
		m = [temp,temp2,temp3,temp4]	
		# print("sum of e = ", (temp+temp2+temp3+n))	
		answer.append(m)
	# print(answer)
	return answer

#make sure its sticky (i.e. >0.9)
def transition():
	answer = []
	for i in range(2):
		m = []		
		temp = random.uniform(0.9,1.0)
		temp2 = 1 - temp
		if (i == 0):
			m.append(temp)
			m.append(temp2)
		else:
			m.append(temp2)
			m.append(temp)
		answer.append(m)
	return answer

#find the actual probability of each state at each time step
def forward_backward(sequence,states,transition_matrix,emission_matrix):
	L = len(sequence)

	forward = []
	f_prev = {}

	#forward algorithm
	for i, x_i in enumerate(sequence):
		f_curr = {}
		for st in states:
			if i == 0:
				# base case for forward
				prev_f_sum = 0.5
			else:
				prev_f_sum = sum(f_prev[k]*transition_matrix[k][st] for k in states)
			f_curr[st] = emission_matrix[st][x_i] * prev_f_sum

		forward.append(f_curr)
		f_prev = f_curr

	p_fwd = sum(f_curr[k]*1.0 for k in states)

	backward = []
	b_prev = {}

	#backward algorithm
	for i, x_i_plus in enumerate(reversed(sequence[1:]+(None,))):
		b_curr = {}
		for st in states:
			if i == 0:
				#base case for backward part
				b_curr[st] = 1.0
			else:
				b_curr[st] = sum((transition_matrix[st][l] * emission_matrix[l][x_i_plus] * b_prev[l]) for l in states)

		backward.insert(0,b_curr)
		b_prev = b_curr

	p_backward = sum((0.5*emission_matrix[l][sequence[0]] * b_curr[l]) for l in states)

	assert (abs(p_fwd - p_backward) < (10**-150))

	return forward,backward,p_fwd

def baumWelch(sequence, states, forward, backward, transition_matrix, emission_matrix):
	
	#Gamma, probability of being in state i at time t
	#calculated in fwd-bwd algorithm(posterior)
	gamma = initialize_E()
	for i in range(len(sequence)):
		for st in states:
			if (sequence[i] == 'A'):
				gamma[st]['A'] += forward[i][st]*backward[i][st]
			if (sequence[i] == 'C'): 
				gamma[st]['C'] += forward[i][st]*backward[i][st]
			if (sequence[i] == 'T'): 
				gamma[st]['T'] += forward[i][st]*backward[i][st]
			if (sequence[i] == 'G'): 
				gamma[st]['G'] += forward[i][st]*backward[i][st]

	#find probability of being in state i and j at times t and t+1
	epsilon = initialize_A()
	for i in range(len(sequence)-1):
		for k in states:
			for l in states:
				epsilon[k][l] += (forward[i][k]* transition_matrix[k][l] * backward[i+1][l] * emission_matrix[l][sequence[i+1]])

	return epsilon,gamma

def print_matrix(matrix):
	for line in matrix:
		print(line)

def print_dict(matrix):
	ans = []
	for i in matrix:
		temp = []
		for j in matrix[i]:
			temp.append(matrix[i][j])
		ans.append(temp)
	print_matrix(ans)

def initialize_E():
	ans = {
		'k' : {'A' : 0.0, 'C' : 0.0, 'T' : 0.0, 'G' : 0.0 },
		'l' : {'A' : 0.0, 'C' : 0.0, 'T' : 0.0, 'G' : 0.0 }
	}
	return ans

def initialize_A():
	ans = {
		'k' : {'k' : 0, 'l' : 0},
		'l' : {'k' : 0, 'l' : 0}
	}
	return ans

def add(x,y):
	for i in x:
		for j in x[i]:
			x[i][j] += y[i][j]	

def divide(x,y):
	for i in x:
		for j in x[i]:
			x[i][j] /= y

def normalize(X):
	for i in X:
		m = 0
		for j in X[i]:
			m += X[i][j]
		for j in X[i]:
			X[i][j] /= m

#about 30 iterations
def run(records,states,transition_probability,emission_probability,prev_log,count):

	log_likelihood = 0
	p_fwd_total = 0
	A = initialize_A()
	E = initialize_E()

	for sequence in records:
		#Forward-Backward algorithm calculation
		forward,backward,p_fwd = forward_backward(tuple(sequence.seq),states,transition_probability,emission_probability)
	
		#perform the baum-welch calculations for A and E
		A_temp,E_temp = baumWelch(sequence,states,forward,backward,transition_probability,emission_probability)

		# print(count,"- Atemp:",A_temp)
		# print(count,"-Etemp",E_temp)
		# print(count," A Before:",A)
		# print(count,"- Ebefore",E)

		#complete summation portion 
		log_likelihood += math.log(p_fwd)
		p_fwd_total += p_fwd
		add(A,A_temp)
		add(E,E_temp)

		# print(count,"A After:",A)
		# print(count,"E After",E)

	divide(A,p_fwd_total)
	divide(E,p_fwd_total)
	
	# print (A)
	# print(E)

	normalize(A)
	normalize(E)

	# print (A)
	# print(E)

	#calculate log_liklihood
	# print(count,":",log_likelihood)
	# print(abs(log_likelihood), abs(prev_log),abs(abs(log_likelihood) - abs(prev_log)),10**-4)
	
	if (abs(abs(log_likelihood) - abs(prev_log)) < 10**-4):
		return transition_probability,emission_probability,log_likelihood
	elif (count > 150):
		print("FAILED!")
		return transition_probability,emission_probability,log_likelihood
	else:
		return run(records,states,A,E,log_likelihood,count+1)

def main():

	#read file and parse into Bio sequence data structure
	records = parse("assign3.fasta")

	"""
	State definitions
	k = A-T rich
	l = G-C rich
	
	Transition Matrix:
	[[kk 	kl]
	 [lk	ll]]

	Emission Matrix:
	[[kA kC kT kG]
	 [lA lC lT lG]
	"""

	states = ('k','l')

	#find normalized initial random emission matrix
	emission_matrix = emission()
	emission_probability = {
		'k' : {'A' : float(emission_matrix[0][0]), 'C' : float(emission_matrix[0][1]), 'T' : float(emission_matrix[0][2]), 'G' : float(emission_matrix[0][3]) },
		'l' : {'A' : float(emission_matrix[1][0]), 'C' : float(emission_matrix[1][1]), 'T' : float(emission_matrix[1][2]), 'G' : float(emission_matrix[1][3]) }
	}

	#find normalized initial "sticky" matrix
	transition_matrix = transition()
	transition_probability = {
		'k' : {'k' : float(transition_matrix[0][0]), 'l' : float(transition_matrix[0][1])},
		'l' : {'k' : float(transition_matrix[1][0]), 'l' : float(transition_matrix[1][1])}
	}

	#test my matrices
	print("\nInitial Transition Matrix: ")
	print_matrix(transition_matrix)
	print("\nInitial Emission Matrix: ")
	print_matrix(emission_matrix)
	print("")

	#calculate forward, backward and baum-welch algorithms for each sequence 
	transition_probability,emission_probability,log_likelihood = run(records,states,transition_probability,emission_probability,1,1)

	print("Log likelihood")
	print(log_likelihood)
	print("\nNew Transition Matrix: ")
	print_dict(transition_probability)
	print("\nNew Emission Matrix: ")
	print_dict(emission_probability)

	return 0

if __name__ == "__main__":
	sys.exit(main())