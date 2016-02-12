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
	indexes = []
	significant_values = []
	answer = (None, None)
	for r in range(len(p)):
		if p[r] < ((r/len(p))*fdr_value):
			answer = (r, p[r])
	return answer

def main():

	#parse in the data file, list of floats
	input_file_name = "data.txt"
	data = parse(input_file_name)
	# print (len(data))

	#plot a histogram of the data
	n, bins, patches = plt.hist(data,bins=40)
	plt.title('Plot of Mixture of Two Different Normal Distributions')
	plt.xlabel('x-values')
	plt.ylabel('y-values')
	plt.show()

	#plot a qq plot of the data
	probplot(data, plot=plt, dist="norm")
	plt.title('QQ Plot of the Data Against the Standard Normal Distribution')
	plt.xlabel('Index')
	plt.ylabel('Normal Data Quantiles')
	plt.show()

	#calculate the p-values for the data
	abs_values = list(abs(n) for n in data)
	abs_values.sort()
	p = list(2.0 * norm.cdf(n * -1) for n in abs_values)
	# print(p)

	#plt a histogram of the p-values
	n, bins, patches = plt.hist(p,bins=40)
	plt.title('Histogram of the P-Values')
	plt.xlabel('Index')
	plt.ylabel('P - Value')
	plt.show()

	#calculate the false discovery rate & the p-value at the FDR cutoff
	p.sort()
	fdr_cutoff = 0.05
	r, p_value = FDR(p, fdr_cutoff)
	print("The p_value at the FDR cutoff of", fdr_cutoff, "is", p_value,".")
	print("There are", r, "p-values less than the pvalue calculated at the FDR cutoff.")

	#plot the raw data values with FDR cutoff
	x = range(len(p))
	plt.scatter(x[:r], p[:r], c='r', edgecolor='r',label='Below FDR cutoff')
	plt.scatter(x[r+1:], p[r+1:], c='b',edgecolor='b',label='Above FDR cutoff')
	# plt.scatter(list(range(len(p_list1),len(p_list1)+len(p_list2))),p_list1)
	plt.plot([0, len(p)], [0.05,0.05], 'k-', lw=2, label="pvalue below 0.05")
	plt.legend(loc='upper left')
	plt.xlim(0,len(p))
	plt.ylim(-.05,1.05)
	plt.xlabel('R - Value')
	plt.ylabel('P - Value')
	plt.title('False Positives of the P-Values')
	plt.show()

	return 0

if __name__ == "__main__":
	sys.exit(main())
