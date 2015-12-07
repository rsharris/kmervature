#!/usr/bin/env python

from sys    import argv,stderr
from string import maketrans


def main():

	assert (len(argv) in [2,3]), "need kmer/count filename and nothing else"

	kmersFilename = argv[1]

	testLimit = None          # we'll only compute FP rate for this many kmers
	if (len(argv) == 3): testLimit = int(argv[2])

	# read kmers file and create kmer-to-count hash

	kmersF = file(kmersFilename, "rt")

	kmerToCount = {}
	for line in kmersF:
		(kmer,count) = line.split()
		kmerToCount[kmer] = int(count)

	kmersF.close()

	# re-read kmers file and count how often a substitution error would be
	# mistaken as a valid kmer

	kmersF = file(kmersFilename, "rt")

	distinctKmers = numMistakable = 0
	totalCount = 0
	weightFP = 0.0

	lineNum = 0
	for line in kmersF:
		lineNum += 1
		if (testLimit != None) and (lineNum > testLimit): break

		kmer = line.split()[0]
		count = kmerToCount[kmer]
		distinctKmers += 1
		totalCount += count

		possibleErrors = kmerFalsePositives = 0
		for errantKmer in errors(kmer):
			possibleErrors += 1
			errorIsFP = (errantKmer in kmerToCount)
			if (not errorIsFP):
				errantKmer = reverse_complement(errantKmer)
				errorIsFP = (errantKmer in kmerToCount)
			if (errorIsFP):
				kmerFalsePositives += 1

		if (kmerFalsePositives > 0): numMistakable += 1
		weightFP += (kmerFalsePositives / float(possibleErrors)) * count
		#print >>stderr, "%s -> %d %d/%d" % (kmer,count,kmerFalsePositives,possibleErrors)

	kmersF.close()

	# report results

	print "distinct kmers = %d" % distinctKmers
	print "distinct kmers with a FP error = %d (%.3f%%)" \
	    % (numMistakable,100.0*numMistakable/distinctKmers)
	print "weighted false positive rate = %.3f%% (%.3f/%d)" \
	    % (100.0*weightFP/totalCount,weightFP,totalCount)


# errors--
#	generate the kmers resulting from a sequencing error

nucToSubs = {"A":"CGT", "C":"AGT", "G":"ACT", "T":"ACG"}

def errors(kmer):
	kmer = [nuc for nuc in kmer]
	for (ix,nuc) in enumerate(kmer):
		for sub in nucToSubs[nuc]:
			kmer[ix] = sub
			yield "".join(kmer)
		kmer[ix] = nuc


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


if __name__ == "__main__": main()
