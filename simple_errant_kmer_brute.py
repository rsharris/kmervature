#!/usr/bin/env python
"""
Compute what the errant kmer distribution can be expected to look like, by
brute force (considering *all* possible error patterns of a given length).

That models the process used in recoverY, though it leaves out such
complications as determining which kmers are good/bad, and that some errant
kmers are actually present elsewhere in the real data.

The output is the distribution of the score-- the fraction of reads expected to
be observed to have a particular score.
"""

from sys import argv,stderr


def main():

	assert (len(argv) == 4), "need kmer_size read_length error_rate and nothing else"

	kmerSize   = int  (argv[1])
	readLength = int  (argv[2])
	errorRate  = float(argv[3])

	assert (0 < kmerSize <= readLength)
	assert (0 <= errorRate < 1)

	assert (readLength <= 20), "max read length is 20"

	# generate all possible reads as good-bad strings

	scoreToP = {}
	for score in xrange(readLength+2-kmerSize):
		scoreToP[score] = 0.0

	for (read,pOfRead) in generate_all_reads(readLength,errorRate):
		score = count_good_kmers(read,kmerSize)
		scoreToP[score] += pOfRead   # $$$ we could suffer precision loss here
		#print pOfRead,read,score

	print "%s\t%s" % ("score_%d_%d" % (kmerSize,readLength),"prob")
	for score in xrange(readLength+2-kmerSize):
		p = scoreToP[score]
		if (p != 0.0): print "%d\t%.6f" % (score,p)
		else:          print "%d\tNA" %  score


def generate_all_reads(readLength,errorRate):
	numReads = 1 << readLength
	pOfNoErrors = (1-errorRate) ** readLength

	for readBits in xrange(numReads):
		read = ["-"] * readLength
		pOfRead = pOfNoErrors

		bits = readBits
		bitNum = 0
		while (bits != 0):
			if ((bits & 1) != 0):
				read[bitNum] = "x"
				pOfRead *= errorRate / (1-errorRate)
			bits >>= 1
			bitNum += 1

		yield ("".join(read),pOfRead)


def count_good_kmers(read,kmerSize):
	goodKmer = "-" * kmerSize
	score = 0
	for ix in xrange(len(read)+1-kmerSize):
		if (read[ix:ix+kmerSize] == goodKmer): score += 1
	return score


if __name__ == "__main__": main()
