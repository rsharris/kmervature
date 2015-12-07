#!/usr/bin/env python
"""
Explore setting of initial parameters for hap-dip model, by trying random
initial parameter sets in the vicinity of the "right" answer, and see where we
lose convergence.
"""

from sys    import argv,stderr
from math   import sqrt
from random import seed as random_seed,random as unit_random
from kmervature import EnrichedHapDipFitter,params_from_text,params_to_text


def main():
	assert (len(argv) == 3), "need the sampleID and number of trials, and nothing else"
	sampleId = argv[1]
	numTrials = int(argv[2])

	random_seed("acorn")
	explainFailure = False
	path = "kmer_histograms"

	# ask the curve fitter what the default paramters are

	fitter = EnrichedHapDipFitter(path+"/"+sampleId+".mixed.kmer_dist")
	paramNames = fitter.paramNames

	defaultParams = fitter.default_params()
	if (defaultParams == None):
		print "(failed to get default params)"
		if (explainFailure):
			print "... return code ..."
			print hdFitter.retCode
			print "... stdout ..."
			print hdFitter.stdout
			print "... stderr ..."
			print hdFitter.stderr
		assert (False)

	defaultParams = params_to_float(defaultParams)

	# read the "good" parameters (usually produced by explore3_hap_dip)

	fitFilename = path+"/"+sampleId+".mixed.fit"

	f = file(fitFilename,"rt")
	goodParams = params_from_text([line for line in f])
	f.close()

	for name in defaultParams:
		assert (name in goodParams), \
		       "parameter \"%s\" missing from %s" % (name,fitFilename)

	for name in goodParams:
		assert (name in defaultParams), \
		       "extra parameter \"%s\" in %s" % (name,fitFilename)

	goodParams = params_to_float(goodParams)

	print params_to_text(paramNames,goodParams,defaultParams,
	                     prefix="good:",prefix2="dflt:")

	# run the convergence trials

	convergenceCount = 0
	for trialNumber in xrange(numTrials):
		print "=== trial %d of %d ===" \
		    % (1+trialNumber,numTrials)

		# choose initial params as a random point in hypercube between "good"
		# and "bad"

		initParams = dict(goodParams)
		norm2Init = 0.0
		for (paramIx,name) in enumerate(paramNames):
			step = unit_random()
			initParams[name] += step*(defaultParams[name]-goodParams[name])
			norm2Init += step*step
		normInit = sqrt(norm2Init) / len(paramNames)

		fitter.set_params(initParams)
		fitParams = fitter.fit()
		if (fitParams == None):
			print params_to_text(paramNames,initParams,prefix="init-[%d]:" % trialNumber)
			print "normInit: %.8f" % normInit
			print "(failure or non-convergence)"
			if (explainFailure):
				print "... return code ..."
				print fitter.retCode
				print "... stdout ..."
				print fitter.stdout
				print "... stderr ..."
				print fitter.stderr
			continue

		print params_to_text(paramNames,initParams,fitParams,
		                     prefix="init+[%d]:" % trialNumber,
		                     prefix2="cvrg[%d]:" % trialNumber)
		fitParams = params_to_float(fitParams)
		dGood = vector_distance(fitParams,goodParams)
		print "normInit: %.8f" % normInit
		print "dGood: %.8f" % dGood
		convergenceCount += 1

	print "%d of %d trials converged" % (convergenceCount,numTrials)


def params_to_float(params):
	return {name:float(params[name]) for name in params}


def vector_distance(vector1,vector2):
	if ([name for name in vector1 if (name not in vector2)] != []): raise valueError
	if ([name for name in vector2 if (name not in vector1)] != []): raise valueError
	return sqrt(sum([((vector1[name]-vector2[name])**2) for name in vector1]))


if __name__ == "__main__": main()
