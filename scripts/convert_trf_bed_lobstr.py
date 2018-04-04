#!/usr/bin/env python

import sys
import argparse

""""
From lobSTR manual:

Column 1: chromosome
Column 2: start coordinate of the STR
Column 3: end coordinate of the STR
Column 4: period of the STR
Column 5: reference copy number
Column 9: STR score. This score measures the purity of the STR sequence and is based on the suggested Tandem Repeats Finder scoring scheme with match=2, mismatch=-7 and indel=-7. Therefore the maximum possible score for a perfectly pure STR sequence (e.g. ATATATATATAT) is 2*(length of STR region).
Column 15: STR repeat unit

TRF manual:

Column 1 and 2: Indices of the repeat relative to the start of the sequence.
Column 3: Period size of the repeat.
Column 4: Number of copies aligned with the consensus pattern.
Column 5: Size of consensus pattern (may differ slightly from the period size).
Column 6: Percent of matches between adjacent copies overall.
Column 7: Percent of indels between adjacent copies overall.
Column 8: Alignment score.
Column 9: Percent composition for each of the four nucleotides.
Column 10: Entropy measure based on percent composition.

The dat file:
Sequence: GmG140828_s1363319



Parameters: 2 7 7 80 10 50 500


930 1002 2 39.5 2 81 15 95 0 0 52 47 1.00 TG TGTGGTGTGGTGGGTGGTGTGTGTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTGTGTGTGTGGTGT
930 1002 11 6.6 12 83 13 95 0 0 52 47 1.00 TGTGTGTGTGTG TGTGGTGTGGTGGGTGGTGTGTGTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTGTGTGTGTGGTGT
2815 2849 2 17.5 2 93 0 61 48 51 0 0 1.00 AC ACACACACACACACACACACACACACACCCACACA


"""


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=
	'Converts a TRF .dat file to a BED suitable for lobSTR.')
	parser.add_argument('-d', '--dat', action='store', help='TRF dat file', required=True, type=argparse.FileType('r'))	

	args = parser.parse_args()
	
	seq = ''
	for line in args.dat:
		split_line = line.strip().split()
		if split_line:
			if split_line[0] == 'Sequence:':
				seq = split_line[1]
			if len(split_line) == 15:
				#consensus pattern can be larger than periode
				if int(split_line[4]) < 7: 
					print ("%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\t%s\tNA\tNA\tNA\tNA\tNA\t%s" % (seq, split_line[0], split_line[1], split_line[2], split_line[3], split_line[7], split_line[13]))
			
			
	

