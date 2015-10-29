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

def main():

	#parse in the PRP file
	filename = "PRP.fasta"
	records = parse(filename)

	print(records)

	return 0

if __name__ == "__main__":
	sys.exit(main())