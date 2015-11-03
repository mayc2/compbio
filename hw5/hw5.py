import sys
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import random
from Bio import SeqIO


def parse(filename):
	records = list(SeqIO.parse(open(filename),"fasta"))
	sequences = []
	ids = []
	for item in records:
		sequences.append(str(item.seq))
		ids.append(item.id)
	return sequences,ids

def calculateConsensus(sequences):
	count_matrix = initConsensus(sequences)
	for i in sequences:
		for j in range(len(i)):
			count_matrix[i[j]][j] += 1.0
	consensus_seq =""
	for j in range(len(sequences[0])):
		consensus_seq = consensus_seq + str(max((count_matrix[key][j],key) for key in ['A','C','U','G','-'])[1])
	return count_matrix, consensus_seq

def initConsensus(sequences):
	answer = {} 
	for key in ['A','C','U','G','-']:
		temp = []
		for _ in range(len(sequences)):
			temp.append(0.0)
		answer[key] = temp
	return answer

def calculateMutualInfo(count_matrix):
	M = initM(count_matrix)

	return M

def initM(count_matrix):
	answer = []
	for _ in range(len(count_matrix['A'])):
		temp = []
		for _ in range(len(count_matrix['A'])):
			temp.append(0.0)
		answer.append(temp)
	return answer

def main():

	#parse in the PRP file
	filename = "alignment.clw"
	sequences,ids = parse(filename)
	# print (sequences)

	
	"""
	Calculate Consensus Sequence
	Matrix dimensions: A,C,T,G,- by position (5xL)
	consensus_seq is a string of max count at each column size
	"""
	count_matrix, consensus_seq = calculateConsensus(sequences)
	# print (count_matrix)
	# print(consensus_seq)

	"""
	Calculate Mutual Information
	Matrix Dimensions: column size x column size
	"""
	M = calculateMutualInfo(count_matrix)
	print(M)

	return 0

if __name__ == "__main__":
	sys.exit(main())