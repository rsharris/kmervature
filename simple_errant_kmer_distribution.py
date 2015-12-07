#!/usr/bin/env python
"""
Simulate what the errant kmer distribution can be expected to look like.

We generate "reads" as simply a sequence of '-" and "x". The latter represents
substituion errors when the read was sequenced.  Then we count how many kmers
do not contain any errors-- the "score" of the read.

That models the process used in recoverY, though it leaves out such
complications as determining which kmers are good/bad, and that some errant
kmers are actually present elsewhere in the real data.

The output is the distribution of the score-- how many reads were observed to
have a particular score.
"""

from sys    import argv,stderr
from random import random as unit_random


def main():

	assert (len(argv) == 5), "need kmer_size read_length error_rate num_reads and nothing else"

	kmerSize   = int  (argv[1])
	readLength = int  (argv[2])
	errorRate  = float(argv[3])
	numReads   = int  (argv[4])

	assert (0 < kmerSize <= readLength)
	assert (0 <= errorRate < 1)

	# run random trials to simulate error distribution

	scoreToCount = {}
	for score in xrange(readLength+2-kmerSize):
		scoreToCount[score] = 0

	for _ in xrange(numReads):
		read  = generate_read(readLength,errorRate)
		score = count_good_kmers(read,kmerSize)
		scoreToCount[score] += 1
		#print >>stderr, "%d %s" % (score,read)

	print "%s\t%s" % ("score_%d_%d" % (kmerSize,readLength),"count")
	for score in xrange(readLength+2-kmerSize):
		count = scoreToCount[score]
		if (count != 0): print "%d\t%d" % (score,count)
		else:            print "%d\tNA" %  score


def generate_read(readLength,errorRate):
	read = ["-"] * readLength
	for ix in xrange(readLength):
		if (unit_random() < errorRate):
			read[ix] = "x"
	return "".join(read)


def count_good_kmers(read,kmerSize):
	goodKmer = "-" * kmerSize
	score = 0
	for ix in xrange(len(read)+1-kmerSize):
		if (read[ix:ix+kmerSize] == goodKmer): score += 1
	return score


if __name__ == "__main__": main()
