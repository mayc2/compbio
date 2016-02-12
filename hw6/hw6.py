from Bio import Phylo
import numpy as np
import matplotlib as plt
import sys
import math

#P(data|tree)
#likelihood of the data given the tree
def felsenstein(tree):
	
	def init_dict():
		answer = {}
		for n in ["A","C","T","G"]:
			answer[n] = 0.0
		return answer

	def init_Sub():
		answer = {}
		for n in ["A","C","T","G"]:
			temp = {}
			for n2 in ["A","C","T","G"]:
				temp[n2] = 0.0
			answer[n] = temp
		return answer

	#4x4 Substitution Matrix
	def jukes_cantor(t):
		alpha = 1.0 * (10 ** (-9))
		t = t * (1.0 * 10 ** 6)
		r_t = (1/4) * (1.0 + 3.0 * math.e ** (-4.0 * alpha * t)) 
		s_t = (1/4) * (1.0 - math.e ** (-4 * alpha * t))
		answer = init_Sub()
		for i in ["A","C","T","G"]:
			for j in ["A","C","T","G"]:
				if i==j:
					answer[i][j] = r_t
				else:
					answer[i][j] = s_t
		# import pprint
		# pprint.pprint(answer, width=1)
		return answer

	#found leaf, base case
	if str(tree) != "Clade":
		P = init_dict()
		P[tree.name[len(tree.name)-1]] = 1.0
		print ("Probability at leaf",tree.name,P)
		return P

	#if not leaf, find children 
	[left,right] = tree.clades	

	D_left = felsenstein(left)
	D_right = felsenstein(right)

	S_left = jukes_cantor(left.branch_length)
	S_right = jukes_cantor(right.branch_length)

	P_left = init_dict()
	P_right = init_dict()
	for n in ["A","C","T","G"]:
		P_left[n] = D_left["A"]*S_left["A"][n]\
			+ D_left["C"]*S_left["C"][n]\
			+ D_left["T"]*S_left["T"][n] \
			+ D_left["G"]*S_left["G"][n]
		P_right[n] = D_right["A"]*S_right["A"][n]\
			+ D_right["C"]*S_right["C"][n]\
			+ D_right["T"]*S_right["T"][n] \
			+ D_right["G"]*S_right["G"][n]

	P_k = init_dict()
	for n in ["A","C","T","G"]:
		P_k[n] = P_left[n]*P_right[n]

	print("Probability at node None", P_k)
	return P_k

def find_likelihood(X):
	answer = 0.0
	for n in ["A","C","T","G"]:
		answer += X[n] * 0.25
	return answer

def parse(filename,filetype):
	return Phylo.read(filename,filetype)

def main():

	#parse the tree files
	filename1 = "tree1.txt"
	filename2 = "tree2.txt"
	filename3 = "tree3.txt"
	filename4 = "tree4.txt"
	filetype = "newick"
	tree1 = parse(filename1,filetype)
	tree2 = parse(filename2, filetype)
	tree3 = parse(filename3, filetype)
	tree4 = parse(filename4, filetype)
	
	#make all trees rooted
	tree1.rooted = True
	tree2.rooted = True
	tree3.rooted = True
	tree4.rooted = True


	#print trees
	print("\nTree 1:\n",tree1,"\n")
	#run recursive fels step
	X1 = felsenstein(tree1.clade)
	#compute likelihood of each calculated probability
	likelihood1 = find_likelihood(X1)
	#print likelihoods
	print("\ntree1's likelihood",likelihood1,"\n")

	print("\nTree 2:\n",tree2,"\n")
	X2 = felsenstein(tree2.clade)
	likelihood2 = find_likelihood(X2)
	print("\ntree2's likelihood",likelihood2,"\n")

	print("\nTree 3:\n",tree3,"\n")
	X3 = felsenstein(tree3.clade)
	likelihood3 = find_likelihood(X3)
	print("\ntree3's likelihood",likelihood3,"\n")

	print("\nTree 4:\n",tree4,"\n")
	X4 = felsenstein(tree4.clade)
	likelihood4 = find_likelihood(X4)
	print("\ntree4's likelihood",likelihood4,"\n")


	return 0

if __name__ == "__main__":
	sys.exit(main())