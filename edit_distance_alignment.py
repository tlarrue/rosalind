'''
CMSC701: Computational Genomics
rosalind problem: reconstructing edit distance

Tara Larrue
'''

import sys
import numpy as np


def read_fasta(input_file):
	'''converts protein strings in FASTA format into list'''

	#read input file
	f = open(input_file, 'rb')
	data = f.read().strip()
	f.close()

	#break up strings into lists		
	pairs = [i.strip(' \n').split('\n') for i in data.split('>')[1:]]

	#construct dictionary
	#dataset = {}
	#for p in pairs:
	#	dataset[p[0]] = p[1]

	dataset = [p[1] for p in pairs]

	return dataset

def compare_symbols(sym1, sym2):
	if sym1 == sym2:
		return 0
	else:
		return 1

def calc_hamming_distance(seq1, seq2):
	'''calculates hamming distance - length of min # of symbols 
	substituions requred to transform one string to another - 
	of 2 protein sequences'''

	#set sequence length & 'gap score'
	m = len(seq1) + 1
	n = len(seq2) + 1
	g = 1;

	#initialize path graph & trackback edge matrix
	hamming = np.zeros((m,n), dtype='i2')
	#traceback_edges = np.zeros((m,n,3), dtype='bool')

	for i in range(m): hamming[i,0] = i*g
	for j in range(n): hamming[0,j] = j*g

	#fill in path graph
	for i in range(1,m):
		for j in range(1,n):

			score_choices = [hamming[i-1,j-1] + compare_symbols(seq1[i-1],seq2[j-1]), #ind0 = match (diag)
							 hamming[i-1,j]+g, #ind1 = gap in seq2 (down arrow)
							 hamming[i,j-1]+g] #ind2 = gap in seq1 (right arrow)

			best_score = min(score_choices)
			
			hamming[i,j] = best_score

			#store best paths
			#is_min = [s==best_score for s in score_choices]
			#traceback_edges[i,j] = is_min

	return hamming[m-1,n-1], hamming


def find_optimal_alignment(hamming_matrix, seq1, seq2):
	'''traces back edges in hamming matrix to find an optimal alignment'''

	optimal1 = [] 
	optimal2 = []

	row = len(seq1)
	col = len(seq2)

	while (row > 0) or (col > 0):

		diag = hamming_matrix[row-1, col-1]
		left = hamming_matrix[row, col-1]
		up = hamming_matrix[row-1,col]

		choices = [diag,up,left]

		ind = choices.index(min(choices))

		if ind == 0:
			optimal1.insert(0,seq1[row-1])
			optimal2.insert(0,seq2[col-1])

			row = row-1
			col = col-1

		elif ind == 1:
			optimal1.insert(0, seq1[row-1])
			optimal2.insert(0, "-")

			row = row - 1

		elif ind == 2:
			optimal1.insert(0, "-")
			optimal2.insert(0, seq2[col-1])

			col = col - 1


	return ''.join(optimal1), ''.join(optimal2)


def main(filename):

	sequences = read_fasta(filename)

	d_h, h_matrix = calc_hamming_distance(*sequences)

	#opt1, opt2 = find_optimal_alignment(t_e, *sequences)
	opt1, opt2 = find_optimal_alignment(h_matrix,*sequences)
	print d_h
	print opt1
	print opt2

if __name__=='__main__':
	input_file = sys.argv[1]
	main(input_file)