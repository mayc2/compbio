import numpy 
from scipy.stats import probplot, norm
import matplotlib.pyplot as plt
from matplotlib import lines
import sys
import math

def parse(input_file_name):
	with open(input_file_name) as f:
		lines = f.read().splitlines()
	for n in range(len(lines)):
		lines[n] = float(lines[n])
	return lines

def FDR(p, fdr_value):
	#find significant p values
	answer = (None, None)
	for r in range(len(p)):
		if p[r][1] < ((r/len(p))*fdr_value):
			answer = (r, p[r][1], p[r][0])
	return answer

def index(data):
	answer = []
	for n in range(len(data)):
		answer.append((n, data[n]))
	return answer

def clean(p):
	answer = []
	for item in p:
		answer.append(item[1])
	return answer

# def sort_p(p):
	# for i in range(len(p)):
		# for 

def main():

	#parse in the data file, list of floats
	input_file_name = "data.txt"
	data = parse(input_file_name)
	# print (len(data), data)

	indexed_data = index(data)
	print(indexed_data)

	#plot a histogram of the data
	n, bins, patches = plt.hist(data,bins=40)
	plt.show()

	#plot a qq plot of the data
	probplot(data, plot=plt, dist="norm")
	plt.show()

	#calculate the p-values for the data
	abs_values = list( (n[0], abs(n[1])) for n in indexed_data)
	abs_values = sort_p(abs_values)
	p = list( (n[0] , 2.0 * norm.cdf(n[1] * -1) ) for n in abs_values)
	print(p)

	p_cleaned = clean(p)

	#plt a histogram of the p-values
	n, bins, patches = plt.hist(p_cleaned,bins=40)
	plt.show()

	#calculate the false discovery rate & the p-value at the FDR cutoff
	p = sort_p(p)
	fdr_cutoff = 0.05
	r, p_value, p_index = FDR(p, fdr_cutoff)
	print(r, p_value)

	#plot the raw data values with FDR cutoff
	plt.plot(list(range(len(data))),data, marker='o', ls='')
	plt.plot([0,len(data)],[p_value,p_value],'k-',lw=2)
	plt.show()

	return 0

if __name__ == "__main__":
	sys.exit(main())
