#!/usr/bin/env python
"""
Compute what the errant kmer distribution can be expected to look like, using a
recurrence formula (formula below).

That models the process used in recoverY, though it leaves out such
complications as determining which kmers are good/bad, and that some errant
kmers are actually present elsewhere in the real data.

The output is the distribution of the score-- the fraction of reads expected to
be observed to have a particular score.

=== Recurrence ===

given:
	k = kmer length
	L = read length
	p = per-base error rate (errors assumed to be idependent)

computes:
	P(S) = probability that a read will have exactly S good kmers

We assume there are no doppleganger kmers (kmers which contain errors but
which exist elsewhere in the sequenced subject).

We define
	P(L,S,r) = probability, among all reads of length L, of having score S with
	           the leftmost error at position r;  r starts at 1;  r>=k
	           indicate the leftmost error is far enough that prepending a non-
	           error will increase S, so we represent any r>=k by r=k.

Conceptually we begin with a read of length 0 with an error immediately
following it.  Though no such position exists, putting an error there keeps us
from counting any kmer as good until the length reaches k.

	P(0,0,1) = 1
	P(0,0,r) = 0 otherwise.

Now extend the length of a general read by 1, on the left.

	P(L,S,r) with probabiliy  p  becomes P(L+1,S,1)
	         with probabiliy 1-p becomes P(L+1,S,r+1) if (r<k) 
	                                or   P(L+1,S+1,k) otherwise

P(S) is just the sum of P(L,S,r) over all r, for L = the read length.
"""

from sys import argv,stderr


def main():

	assert (len(argv) == 4), "need kmer_size read_length error_rate and nothing else"

	kmerSize   = int  (argv[1])
	readLength = int  (argv[2])
	errorRate  = float(argv[3])

	assert (0 < kmerSize <= readLength)
	assert (0 <= errorRate < 1)

	# compute the recurrence

	k = kmerSize
	p = errorRate
	q = 1-p

	P = {}

	L = 0
	P[L] = {}
	P[L][(0,1)] = 1.0

	for L in xrange(0,readLength):
		P[L+1] = {}

		for (S,r) in P[L]:
			pShort = P[L][(S,r)]

			# prepend error

			outcome = (S,1)
			if (outcome not in P[L+1]): P[L+1][outcome] = 0.0
			P[L+1][outcome] += p*pShort

			# prepend non-error

			if (r < k):
				outcome = (S,r+1)
				if (outcome not in P[L+1]): P[L+1][outcome] = 0.0
				P[L+1][outcome] += q*pShort
			else:
				outcome = (S+1,k)
				if (outcome not in P[L+1]): P[L+1][outcome] = 0.0
				P[L+1][outcome] += q*pShort

	# compute (by summing) and print the distribution for the final length
	# nota bene: we could print the distribution for every length

	L = readLength
	PSUM = {}
	for (S,r) in P[L]:
		if (S not in PSUM): PSUM[S] = 0.0
		PSUM[S] += P[L][(S,r)]

	print "%s\t%s" % ("score_%d_%d" % (kmerSize,readLength),"prob")
	for S in xrange(readLength+2-kmerSize):
		if (S not in PSUM): print "%d\tNA" %  S
		else:               print "%d\t%.6f" % (S,PSUM[S])


if __name__ == "__main__": main()
