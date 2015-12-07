#!/usr/bin/env python
"""
Compute the ground truth of kmer counts for fake reads sampled from a given
sequence.
"""

from sys    import argv,stdin,stdout,stderr,exit
from string import maketrans
try:                from hashlib import md5 as md5_new
except ImportError: from md5     import new as md5_new


def usage(s=None):
	message = """
usage: cat reads | kmer_ground_truth <kmer_size> <source_fasta> [options]
  <kmer_size>           (required) number of nucleotides in a kmer
  <source_fasta>        (required) the sequence the reads were sampled from
  <source_fasta2>       (optional) if the reads were sampled from a diploid
                        sequence, this is that sequence; there are certain
                        expectations of sequence name format in this case (see
                        below)
  M=[<residue>/]<modulus> screen kmers with a residue;  first number is the
                        residue (plus 1), second is the modulus; this is a
                        memory-reduction technique; if the residue is absent,
                        1 is used by default
  --kmerlimit=none      don't limit the number of "distinct" kmers we'll handle
  --kmerlimit=<number>  limit the number of "distinct" kmers we'll handle; if
                        this limit is exceded, we quit (and don't write the
                        histograms); this is a guard against using excessive
                        memory and some kmers may be counted more than once
                        (default is 1G).
  --head=<number>       limit the number of input reads we'll process
  --progress=<number>   periodically report how many reads we've processed

Compute the ground truth of kmer counts for fake reads sampled from a given
sequence or sequences.

We expect reads to have a name that indicates where they were sampled from, in
the form [{anyprefix}_]{chrom}_{zstart}_{strand}.
	{anyprefix} is option and can be anything
	{chrom}     is typically the name of the chromosome; beware that this can't
	            contain an underscore; for diploid processing, we assume
	            {chrom} is of the form {name}.fragment{1|2}
	{zstart}    is the zero-based start of the interval the read was sampled
	            from, always counted along the forward strand
	{strand}    is either F (farward) or R (reverse-complement)"""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global debug

	# parse the command line

	kmerSize          = None
	sourceFilename    = None
	source2Filename   = None
	modulus = residue = None
	kmerLimit         = int_with_unit("1G")
	headLimit         = None
	reportProgress    = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("M=")) or (arg.startswith("--subspace=")):
			if (not "/" in argVal):
				modulus = int_with_unit(argVal)
				residue = 1
				assert (0 < modulus)
			else:
				(residue,modulus) = argVal.split("/",1)
				modulus = int_with_unit(modulus)
				residue = int_with_unit(residue)
				assert (0 < residue <= modulus)
			if (modulus == 1) and (residue == 1):
				modulus = residue = None
		elif (arg.startswith("--kmerlimit=")):
			if (argVal.lower() == "none"): kmerLimit = None
			else:                          kmerLimit = int_with_unit(argVal)
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (kmerSize == None):
			try:
				kmerSize = int(arg)
				if (kmerSize < 1): raise ValueError
			except ValueError:
				usage("invalid kmer size: \"%s\"" % arg)
		elif (sourceFilename == None):
			sourceFilename = arg
		elif (source2Filename == None):
			source2Filename = arg
		else:
			usage("unrecognized option: %s" % arg)

	# read the source(s)

	f = file(sourceFilename,"rt")
	source = {}
	for (chrom,seq) in read_fasta(f):
		source[chrom] = seq
	f.close()
	#print [(chrom,len(source[chrom])) for chrom in source]

	source2 = None
	if (source2Filename != None):
		f = file(source2Filename,"rt")
		source2 = {}
		for (chrom,seq) in read_fasta(f):
			source2[chrom] = seq
		f.close()
		#print [(chrom,len(source2[chrom])) for chrom in source2]

	# process the reads, to collect abundance counts of good and bad kmers

	goodKmerToAbundance = {}
	errKmerToAbundance  = {}
	hetKmerToAbundance  = {}
	numDistinctKmers = 0

	goodKmerSyndrome = "-" * kmerSize

	readNumber = 0
	for (name,read) in read_fasta(stdin):
		readNumber += 1

		if (reportProgress != None) and (readNumber % reportProgress == 1) and (readNumber != 1):
			progressCount = commatize(readNumber-1)
			if (kmerLimit == None):
				print >>stderr, "progress: %s reads processed" % progressCount
			else:
				print >>stderr, "progress: %s reads processed (%s \"distinct\" %d-mers)" \
			     % (progressCount,commatize(numDistinctKmers),kmerSize)

		if (headLimit != None) and (readNumber > headLimit):
			print >>stderr, "limit of %s reads reached" % commatize(headLimit)
			readNumber -= 1
			break

		assert (kmerLimit == None) or (numDistinctKmers <= kmerLimit), \
		       "limit of %s \"distinct\" %d-mers exceeded (%s reads)" \
		     % (commatize(kmerLimit),kmerSize,commatize(readNumber-1))

		(chrom,start,strand) = parse_read_name(name)

		if (chrom in source):
			seq = source[chrom]
			srcIs2 = False
		elif (source2 != None) and (chrom in source2):
			seq = source2[chrom]
			srcIs2 = True
		else:
			assert (False), "no source was provided for \"%s\" (read \"%s\")" \
			              % (chrom,name)

		end = start + len(read)
		if (end > len(seq)):
			assert (False), "%d..%d is beyond end of \"%s\" (read \"%s\")" \
			              % (start,end,chrom,name)

		src = seq[start:end]
		if (strand == "R"): src = reverse_complement(src)
		syndrome = reduce_to_mismatches(read,src)

		if (source2 != None):
			if (srcIs2):
				chromOther = chrom.replace("fragment2","fragment1")
				assert (chromOther in source), \
				       "no source was provided for \"%s\" (parallel to read \"%s\")" \
					 % (chromOther,name)
				seqOther = source[chromOther]
			else:
				chromOther = chrom.replace("fragment1","fragment2")
				assert (chromOther in source2), \
				       "no source was provided for \"%s\" (parallel to read \"%s\")" \
					 % (chromOther,name)
				seqOther = source2[chromOther]
			srcOther = seqOther[start:end]
			if (strand == "R"): srcOther = reverse_complement(srcOther)

		for ix in xrange(len(syndrome)+1-kmerSize):
			kmer  = read[ix:ix+kmerSize]
			canon = canonical_kmer(kmer)

			if (modulus != None):
				if (hash_of_kmer(canon,modulus) != residue): continue

			isHeterozygous = False
			if (source2 != None):
				kmerOther = srcOther[ix:ix+kmerSize]
				isHeterozygous = (kmerOther != kmer)

			if (syndrome[ix:ix+kmerSize] != goodKmerSyndrome):
				if (canon not in errKmerToAbundance):
					errKmerToAbundance[canon] =  1
					numDistinctKmers += 1
				else:
					errKmerToAbundance[canon] += 1
			else:
				if (canon not in goodKmerToAbundance):
					goodKmerToAbundance[canon] =  1
					numDistinctKmers += 1
				else:
					goodKmerToAbundance[canon] += 1
				if (isHeterozygous):
					if (canon not in hetKmerToAbundance):
						hetKmerToAbundance[canon] =  1
						numDistinctKmers += 1
					else:
						hetKmerToAbundance[canon] += 1

	# report results

	goodAbundanceToCount = {}
	for kmer in goodKmerToAbundance:
		a = goodKmerToAbundance[kmer]
		if (a not in goodAbundanceToCount): goodAbundanceToCount[a] =  1
		else:                               goodAbundanceToCount[a] += 1

	errAbundanceToCount = {}
	for kmer in errKmerToAbundance:
		a = errKmerToAbundance[kmer]
		if (a not in errAbundanceToCount): errAbundanceToCount[a] =  1
		else:                              errAbundanceToCount[a] += 1

	abundances =  [a for a in goodAbundanceToCount]
	abundances += [a for a in errAbundanceToCount if (a not in abundances)]
	abundances.sort()

	if (source2 != None):
		hetAbundanceToCount = {}
		homAbundanceToCount = {}
		for kmer in goodKmerToAbundance:
			if (kmer in hetKmerToAbundance): a = hetKmerToAbundance[kmer]
			else:                            a = 0
			if (a not in hetAbundanceToCount): hetAbundanceToCount[a] =  1
			else:                              hetAbundanceToCount[a] += 1
			a = goodKmerToAbundance[kmer] - a
			if (a not in homAbundanceToCount): homAbundanceToCount[a] =  1
			else:                              homAbundanceToCount[a] += 1

	if (source2 == None):
		print "%s\t%s\t%s" % ("abundance","good","error")
	else:
		print "%s\t%s\t%s\t%s\t%s" % ("abundance","good","error","homozygous","heterozygous")

	for a in abundances:
		g = b = "NA"
		if (a in goodAbundanceToCount): g = goodAbundanceToCount[a]
		if (a in errAbundanceToCount):  b = errAbundanceToCount[a]

		if (source2 == None):
			print "%d\t%s\t%s" % (a,g,b)
			continue

		hom = het = "NA"
		if (a in hetAbundanceToCount): het = hetAbundanceToCount[a]
		if (a in homAbundanceToCount): hom = homAbundanceToCount[a]

		print "%d\t%s\t%s\t%s\t%s" % (a,g,b,hom,het)

	if (kmerLimit != None):
		print >>stderr, "%s reads had %s \"distinct\" %d-mers" \
	     % (commatize(readNumber),commatize(numDistinctKmers),kmerSize)


# reduce_to_mismatches--
#	compute the mismatch string of a read and the source it was drawn from

def reduce_to_mismatches(read,src):
	readLength = len(read)
	syndrome = ["-"] * readLength
	for (ix,readNuc) in enumerate(read):
		if (readNuc != src[ix]): syndrome[ix] = "x"
	return "".join(syndrome)


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


# parse_read_name--

def parse_read_name(name):
	tokens = name.split("_")
	try:
		if (len(tokens) < 3): raise ValueError
		(chrom,start,strand) = tokens[-3:]
		start = int(start)
		if (start < 0): raise ValueError
		if (strand not in ["F","R"]): raise ValueError
	except ValueError:
		assert (False), "can't parse read name \"%s\"" % name
	return (chrom,start,strand)


# canonical_kmer--
#	Choose a consistent representative for any kmer and its reverse compliment.

def canonical_kmer(kmer):
	rev = reverse_complement(kmer)
	if (kmer < rev): return kmer
	else:            return rev


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# hash_of_kmer--
#	$$$ md5 is overkill for this, and we'd like to use something faster, but
#	$$$ .. str.__hash__() isn't very good for this purpose
#	Reduce a kmer to a hash value modulo some modulus.  The value h returned
#	is in the range 0 < h <= modulus

def hash_of_kmer(kmer,modulus):
	h = md5_new()
	h.update(kmer)
	return 1 + (int(h.hexdigest()[:25],16) % modulus)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


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
