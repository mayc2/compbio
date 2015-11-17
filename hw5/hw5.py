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

def convert(count_matrix):
	L = len(count_matrix['A'])

	frequency_matrix = count_matrix
	for key in ['A','C','U','G','-']:
		for i in range(len(frequency_matrix[key])):
			frequency_matrix[key][i] /= L
	return frequency_matrix

def calculateMutualInfo(frequency_matrix,sequences):

	#initialize matrix
	L = len(frequency_matrix['A'])
	M = np.zeros([L,L])

	def combined_freq(i,j,n1,n2):
		answer = 0.0
		for seq in sequences:
			if seq[i] == n1 and seq[j] == n2:
				answer += 1
		return answer

	def find1(i,j):
		answer = 0.0
		for n1 in ['A','C','G','U','-']:
			for n2 in ['A','C','G','U','-']:
				if delta(n1,n2) == 1:
					answer += combined_freq(i,j,n1,n2)
		return answer

	def find2(i,j):
		answer = 0.0
		for n1 in ['A','C','G','U','-']:
			for n2 in ['A','C','G','U','-']:
				f_i = frequency_matrix[n1][i]
				f_j = frequency_matrix[n2][j]
				f_ij = 0.0
				if delta(n1,n2) == 1:
					f_ij = combined_freq(i,j,n1,n2)
				if (f_i != 0 and f_j != 0 and f_ij != 0):	
					answer += math.log(f_ij / (f_i *f_j))
		return answer 

	for i in range(L):
		for j in range(L):
			if i != j:
				value1 = find1(i,j)
				value2 = find2(i,j)
				answer = value1 * value2	
				if answer < 0:
					M[i][j] = 0.0
				else:
					M[i][j] = answer
	return M

def delta(si, sj):
    """
    delta - score base pairs for Nussinov algorithm
    si, si - nucleotide letters
    return
    d - delta score 1 for complement or G-U wobble pair, otherwise 0
        1 indicates possible base pairng
    """
    
    nt2int = {'A': 0, 'C': 1, 'G': 2, '-': 4, 'U': 3}
    nt2comp = {'A': 3, 'C': 2, 'G': 1, '-': 5, 'U': 0}

    d = 0
    if nt2int[si] == nt2comp[sj]:
        d = 1
    elif (si == 'G' and sj == 'U') or (si == 'U' and sj == 'G'):
        d = 1

    return d

def dynamic_program(M):
	
	#initiaize matrix
	L = len(M)
	D = np.zeros([L,L])

	def find4(i,j):
		answer = 0.0
		k = i + 4
		while(k <= (j-5)): 
			answer = max(answer, D[i][k]+D[k+1][j])
			# print(answer)
			k += 1
		return answer

	#recursion
	for i, e in reversed(list(enumerate(M))):
		for j in range(i+3,L):
			value1 = D[i+1][j]
			value2 = D[i][j-1]
			value3 = D[i+1][j-1]+M[i][j]
			value4 = find4(i,j)
			# print(value1,value2,value3,value4)
			D[i][j] = max(value1,value2,value3,value4)
	return D

def BackTrace(g, M):
    """
    BackTrace - backtrace throug Nussinov matrix
    g - matrix returned by Nussinov
    seq - sequence string
    return
    pairs - pairs array; i, i+1 base paired
    M - numpy pairs array M[i,j] == 1 if i and j base paired
    """

    def traceback(i, j):
        """
        traceback - recursive back trace through g matrix
                    fills pairs matrix defined in outer scope
        """
        nonlocal pairs
        
        if i < j:
            if g[i,j] == g[i+1, j]:
                traceback(i+1, j)
            elif g[i,j] == g[i, j-1]:
                traceback(i, j-1)
            elif g[i,j] == g[i+1, j-1] + M[i][j]:
                if abs(i-j) > 2:
                	pairs += [i,j]
                traceback(i+1, j-1)
            else:
                for k in range(i+1, j):
                    if g[i,j] == g[i, k] + g[k+1, j]:
                        traceback(i, k)
                        traceback(k+1, j)
                        break
    
    L = g.shape[0]
    pairs = []
    traceback(0, L-1)

    dot_bracket = np.zeros([L, L], dtype =int)
    for i in range(0, len(pairs), 2):
        dot_bracket[pairs[i], pairs[i+1]] = 1

    return pairs, dot_bracket

def StructureFromPairs(pairs, L):
    struct = list('.' * L)
    for i in range(0, len(pairs), 2):
        struct[pairs[i]] = '('
        struct[pairs[i+1]] = ')'
        
    return ''.join(struct)

def main():
	np.set_printoptions(threshold=np.nan)
	#parse in the PRP file or alignment file
	filename = "alignment.clw"
	# filename = "PRP.fasta"
	sequences,ids = parse(filename)
	# print (sequences)

	
	"""
	Calculate Consensus Sequence
	Matrix dimensions: A,C,T,G,- by position (5xL)
	consensus_seq is a string of max count at each column size
	"""
	count_matrix, consensus_seq = calculateConsensus(sequences)
	# print (count_matrix)
	print("The consensus sequence is:")
	print(consensus_seq)

	"""
	Convert count_matrix to a frequency_matrix
	"""
	frequency_matrix = convert(count_matrix)
	# print(frequency_matrix)

	"""
	Calculate Mutual Information
	Matrix Dimensions: column size x column size
	"""
	M = calculateMutualInfo(frequency_matrix,sequences)
	# print(M)

	"""
	Calculate Maximum mutual information secondary structure using dynamic programming
	"""
	g = dynamic_program(M)
	# print(g)

	"""
	Run back trace algorithm to find pairs and the secondary structure in dot-bracket format
	"""
	pairs, dot_bracket = BackTrace(g,M)
	# print(pairs)
	struct = StructureFromPairs(pairs, len(M))
	print("\nThe dot-bracket format description of the consensuse sequence is:")
	print(struct)

	"""
	Plot the 2d secondary structure w/ external source
	"""

	return 0

if __name__ == "__main__":
	sys.exit(main())