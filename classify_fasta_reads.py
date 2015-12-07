#!/usr/bin/env python

from sys    import argv,stderr,getsizeof
from string import maketrans

debug = []
debug += ["memory"]
#debug += ["profile"]
#debug += ["readscore"]


def main():
	global goodKmers,kmerSize

	assert (len(argv) == 4), \
	       "need kmer filename, reads filename, chromosome name, and nothing else"

	kmersFilename    = argv[1]
	fastaFilename    = argv[2]
	targetChromosome = argv[3]

	# read kmers file and create good kmers set
	# $$$ we assume all kmers are the same length, but we don't verify that

	kmersF = file(kmersFilename, "rt")
	goodKmers = set(line.split()[0] for line in kmersF)
	kmersF.close()

	assert (len(goodKmers) != 0)

	for kmer in goodKmers:
		kmerSize = len(kmer)
		break					# (only bother with the first kmer)

	if ("memory" in debug):
		for kmer in goodKmers:
			kmerBytes = getsizeof(kmer)
			break				# (only bother with the first kmer)

		print >>stderr, "each %d-mer uses %s bytes" \
		              % (kmerSize,kmerBytes)
		print >>stderr, "%s %d-mers is %s bytes" \
		              % (commatize(len(goodKmers)),kmerSize,
		                 commatize(len(goodKmers)*kmerBytes))
		print >>stderr, "goodKmers uses %s bytes" \
		              % commatize(getsizeof(goodKmers))

	# process the reads
	# $$$ we assume all reads are the same length, but we don't verify that

	scoreToCount          = {}
	scoreToCountTarget    = {}
	scoreToCountNonTarget = {}

	fastaF = file(fastaFilename, "rt")
	for (name,seq) in read_fasta(fastaF):
		score = count_good_kmers(seq)
		if (score not in scoreToCount):
			scoreToCount[score]          = 0
			scoreToCountTarget[score]    = 0
			scoreToCountNonTarget[score] = 0
		scoreToCount[score] += 1

		if (targetChromosome in name): scoreToCountTarget[score] += 1
		else:                          scoreToCountNonTarget[score] += 1

		if ("readscore" in debug):
			print >>stderr, "%s %d %s" % (targetChromosome in name,score,name)

	fastaF.close()

	# report results

	print "%s\t%s\t%s\t%s" % ("score","observed","target","nontarget")

	scores = [score for score in scoreToCount]
	scores.sort()
	for score in scores:
		print "%d\t%d\t%d\t%d" \
		    % (score,scoreToCount[score],scoreToCountTarget[score],scoreToCountNonTarget[score])


# count_good_kmers--

def count_good_kmers(seq):
	reverseSeq = reverse_complement(seq)

	forwardIx = 0
	reverseIx = len(seq) - kmerSize
	numGood = 0
	profile = []

	while (forwardIx <= (len(seq)-kmerSize)):
		forwardKmer = seq[forwardIx:forwardIx+kmerSize]
		reverseKmer = reverseSeq[reverseIx:reverseIx+kmerSize]
		if ((forwardKmer in goodKmers) or (reverseKmer in goodKmers)):
			numGood += 1
			profile += ["-"]
		else:
			profile += ["x"]
		forwardIx += 1
		reverseIx -= 1

	if ("profile" in debug):
		print >>stderr, name
		print >>stderr, seq
		print >>stderr, "".join(profile)
		print >>stderr

	return numGood


# read_fasta--

def read_fasta(f):
	name = "(nameless)"
	seq  = []

	for line in f:
		line = line.rstrip()
		if (line.startswith(">")):
			if (seq != []): yield (name,"".join(seq))
			name = line[1:].split()[0]
			seq  = []
		else:
			seq += [line]

	if (seq != []): yield (name,"".join(seq))


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# commatize--
#	Convert a numeric string into one with commas.

def commatize(s):
	if (type(s) != str): s = str(s)
	(prefix,val,suffix) = ("",s,"")
	if (val.startswith("-")): (prefix,val) = ("-",val[1:])
	if ("." in val):          (val,suffix) = val.split(".",1)

	try:    int(val)
	except: return s

	digits = len(val)
	if (digits > 3):
		leader = digits % 3
		chunks = []
		if (leader != 0):
			chunks += [val[:leader]]
		chunks += [val[ix:ix+3] for ix in xrange(leader,digits,3)]
		val = ",".join(chunks)

	return prefix + val + suffix


if __name__ == "__main__": main()
