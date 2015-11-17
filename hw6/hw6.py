import Bio
import numpy as np
import matplotlib as plt
import sys

def main():

	filename1 = "tree1.txt"
	filename2 = "tree2.txt"
	filename3 = "tree3.txt"
	tree1 = Bio.Phylo.read(filename1)
	tree2 = Bio.Phylo.read(filename2)
	tree3 = Bio.Phylo.read(filename3)

	print(tree1)
	print(tree2)
	print(tree3)


	return 0

if __init__ == "__main__":
	sys.exit(main())