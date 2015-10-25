import numpy
import scipy
import math
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
import random
import numpy as np

#parse input file into record data structure
def parse():
	chr12 = SeqIO.read(open("chr12.fasta"), "fasta")
	records = list(SeqIO.parse("CpG","fasta"))
	test_sequences = list(SeqIO.parse("test_sequences.fasta","fasta"))
	return records, chr12, test_sequences

#check to make sure selected sequence does not overlap any previous windows
def check_window(prev_windows,sample_length,limits):
	start = random.randrange(limits[0],limits[1]) 
	end = start + sample_length
	
	for window in prev_windows:
		length = window[1] - window[0]
		if (start >= (window[0] - length) and start <= (window[1] + length)):
			start, end = check_window(prev_windows,sample_length,limits)
		if (end >= (window[0] - length) and end <= (window[1] + length)):
			start, end = check_window(prev_windows,sample_length,limits)
		if (end > limits[1]):
			start, end = check_window(prev_windows,sample_length,limits)
	return start, end	

#randomly select another sequence
def select_seq(samples,records,prev_windows,limits,n):
	answer = []
	for i in range(0,n):

		start, end = check_window(prev_windows, len(samples[i]), limits)
		prev_windows.append((start,end))
		answer.append(records.seq[start:end])

	return answer

#find transition matrix from a list of strings
def transition(record):
	a_length = 0.0
	c_length = 0.0
	t_length = 0.0
	g_length = 0.0
	answer = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
	for seq in record:
		for i in range(1,len(seq)):
			#current is A
			if(seq[i-1] == "A" or seq[i-1] == "a"): 
				a_length += 1
				if(seq[i] == "A" or seq[i] == "a"):
					answer[0][0] += 1
				if(seq[i] == "C" or seq[i] == "c"):
					answer[0][1] +=1
				if(seq[i] == "T" or seq[i] == "t"):
					answer[0][2] +=1
				if(seq[i] == "G" or seq[i] == "g"):
					answer[0][3] +=1
			#current is C
			if(seq[i-1] == "C" or seq[i-1] == "c"): 
				c_length += 1
				if(seq[i] == "A" or seq[i] == "a"):
					answer[1][0] += 1
				if(seq[i] == "C" or seq[i] == "c"):
					answer[1][1] +=1
				if(seq[i] == "T" or seq[i] == "t"):
					answer[1][2] +=1
				if(seq[i] == "G" or seq[i] == "g"):
					answer[1][3] +=1
			#current is T
			if(seq[i-1] == "T" or seq[i-1] == "t"): 
				t_length += 1
				if(seq[i] == "A" or seq[i] == "a"):
					answer[2][0] += 1
				if(seq[i] == "C" or seq[i] == "c"):
					answer[2][1] +=1
				if(seq[i] == "T" or seq[i] == "t"):
					answer[2][2] +=1
				if(seq[i] == "G" or seq[i] == "g"):
					answer[2][3] +=1
			#current is G
			if(seq[i-1] == "G" or seq[i-1] == "g"): 
				g_length += 1
				if(seq[i] == "A" or seq[i] == "a"):
					answer[3][0] += 1
				if(seq[i] == "C" or seq[i] == "c"):
					answer[3][1] +=1
				if(seq[i] == "T" or seq[i] == "t"):
					answer[3][2] +=1
				if(seq[i] == "G" or seq[i] == "g"):
					answer[3][3] +=1
	
	for i in range(0,len(answer[0])):
		answer[0][i] /= a_length
		answer[1][i] /= c_length
		answer[2][i] /= t_length
		answer[3][i] /= g_length
	return answer

#output formatted matrix
def print_matrix(matrix):
	print ("         A           C           T           G")
	for i in range(0,len(matrix)):
		if (i == 0):
			print ("A",matrix[i])
		if (i == 1):
			print("C",matrix[i])
		if (i == 2):
			print("T",matrix[i])
		if (i == 3):
			print("G",matrix[i])

def log_odds_ratio(positive_matrix,negative_matrix,records):
	answer = []
	for seq in records:
		temp = 0.0
		for i in range(1,len(seq)):
			#current is A
			if(seq[i-1] == "A" or seq[i-1] == "a"): 
				if(seq[i] == "A" or seq[i] == "a"):
					temp += math.log(positive_matrix.item((0,0))/negative_matrix.item((0,0)))
				if(seq[i] == "C" or seq[i] == "c"):
					temp += math.log(positive_matrix.item((0,1))/negative_matrix.item((0,1)))
				if(seq[i] == "T" or seq[i] == "t"):
					temp += math.log(positive_matrix.item((0,2))/negative_matrix.item((0,2)))
				if(seq[i] == "G" or seq[i] == "g"):
					temp += math.log(positive_matrix.item((0,3))/negative_matrix.item((0,3)))
			#current is C
			if(seq[i-1] == "C" or seq[i-1] == "c"): 
				if(seq[i] == "A" or seq[i] == "a"):
					temp += math.log(positive_matrix.item((1,0))/negative_matrix.item((1,0)))
				if(seq[i] == "C" or seq[i] == "c"):
					temp += math.log(positive_matrix.item((1,1))/negative_matrix.item((1,1)))
				if(seq[i] == "T" or seq[i] == "t"):
					temp += math.log(positive_matrix.item((1,2))/negative_matrix.item((1,2)))
				if(seq[i] == "G" or seq[i] == "g"):
					temp += math.log(positive_matrix.item((1,3))/negative_matrix.item((1,3)))
			#current is T
			if(seq[i-1] == "T" or seq[i-1] == "t"): 
				if(seq[i] == "A" or seq[i] == "a"):
					temp += math.log(positive_matrix.item((2,0))/negative_matrix.item((2,0)))
				if(seq[i] == "C" or seq[i] == "c"):
					temp += math.log(positive_matrix.item((2,1))/negative_matrix.item((2,1)))
				if(seq[i] == "T" or seq[i] == "t"):
					temp += math.log(positive_matrix.item((2,2))/negative_matrix.item((2,2)))
				if(seq[i] == "G" or seq[i] == "g"):
					temp += math.log(positive_matrix.item((2,3))/negative_matrix.item((2,3)))
			#current is G
			if(seq[i-1] == "G" or seq[i-1] == "g"): 
				if(seq[i] == "A" or seq[i] == "a"):
					temp += math.log(positive_matrix.item((3,0))/negative_matrix.item((3,0)))
				if(seq[i] == "C" or seq[i] == "c"):
					temp += math.log(positive_matrix.item((3,1))/negative_matrix.item((3,1)))
				if(seq[i] == "T" or seq[i] == "t"):
					temp += math.log(positive_matrix.item((3,2))/negative_matrix.item((3,2)))
				if(seq[i] == "G" or seq[i] == "g"):
					temp += math.log(positive_matrix.item((3,3))/negative_matrix.item((3,3)))
		temp /= float(len(seq))
		answer.append(temp)
	return answer

def plot_log_odds(positive_log_odds, negative_log_odds):
	plt.hist(positive_log_odds, bins=20, label="CpG")
	plt.hist(negative_log_odds, bins=20, label="Random Chr12")
	plt.title("hg38 CpG Log-Odds Ratios")
	plt.ylabel("Counts")
	plt.xlabel("Log Odds")
	plt.legend()
	plt.show()

def main():
	records, chr12, temp = parse()

	#converting test seq to list of strings for easier data manipulation
	test_sequences = []
	for seq in temp:
		test_sequences.append(seq.seq)
	del temp

	#defined parameters of window to select from in chr12 
	limits = (10200,133265309)

	#variable to track windows that have previously been used
	prev_windows = []

	#find the random 200 positive cpg test set
	temp = random.sample(records,200)
	#remove 200 cpg to get the remaining 1011 cpg positive test set
	cpg_200 = []
	for sample in temp:
		records.remove(sample)
		cpg_200.append(sample.seq)

	#clearing memory for no longer needed list of 200 seq items
	del temp

	#find the negative 200 test set
	negative_200 = select_seq(cpg_200,chr12,prev_windows,limits,200)
	
	#convert 1011 cpg positive to list of strings
	positive_1011 = []
	for seq in records:
		positive_1011.append(seq.seq)

	#find the remaining 1011 negative test set
	negative_1011 = select_seq(records,chr12,prev_windows,limits,1011)

	#clearing memory for initial records list that is no longer needed
	del records

	#find transtion matrices of 200 test sets
	positive_transition = np.matrix(transition(cpg_200))
	negative_transition = np.matrix(transition(negative_200))

	"""
	Transition Matrix Layout
	[[AA  AC  AT  AG]
	 [CA  CC  CT  CG]
	 [TA  TC  TT  TG]
	 [GA  GC  GT  GG]]

	"""
	print("Positive Transition Matrix")
	print_matrix(positive_transition)
	print("")
	print("Negative Transition Matrix")
	print_matrix(negative_transition)	

	positive_log_odds = log_odds_ratio(positive_transition, negative_transition, positive_1011)
	negative_log_odds = log_odds_ratio(positive_transition, negative_transition, negative_1011)

	plot_log_odds(positive_log_odds,negative_log_odds)
	
	test_log_odds = log_odds_ratio(positive_transition,negative_transition,test_sequences)

	print("\nTest Sequence Log Odds Ratios:")
	for i in range(0,len(test_log_odds)):
		print (i,": ",test_log_odds[i])

	# for i in range(0,5):
	# 	print ""
	# 	print (negative_200[i])

if __name__ == "__main__":
	main()