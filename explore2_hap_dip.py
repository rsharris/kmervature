#!/usr/bin/env python
"""
(This functionality has been replaced by explore4_hap_dip)

Explore setting of initial parameters for hap-dip model, by trying random
initial parameter sets in the vicinity of the "right" answer, and see where we
lose convergence.
"""

from sys    import argv,stderr
from math   import sqrt
from random import seed as random_seed,random as unit_random
from kmervature import EnrichedHapDipFitter,params_to_text


def main():
	assert (len(argv) == 1), "give me no arguments"

	numTrials = 1000
	random_seed("acorn")
	explainFailure = False
	path = "kmer_histograms"

	#sampleId = "mixedB"
	#defaultParams = {"zp.copy.y"   :  3.000,
	#                 "zp.copy.hom" :  3.000,
	#                 "zp.copy.het" :  3.000,
	#                 "p.e"         :  0.942,
	#                 "shape.e"     :  3.000,
	#                 "scale.e"     :  1.000,
	#                 "p.y"         :  0.900,
	#                 "u.y"         : 64.000,
	#                 "sd.y"        : 14.826,
	#                 "shape.y"     :  0.000,
	#                 "p.hom"       :  0.800,
	#                 "u.hom"       :  5.120,
	#                 "sd.hom"      :  1.186,
	#                 "var.het"     :  1.407}
	#goodParams    = {"zp.copy.y"   :  2.042,
	#                 "zp.copy.hom" :  3.157,
	#                 "zp.copy.het" : 17.795,
	#                 "p.e"         :  0.935,
	#                 "shape.e"     :  0.096,
	#                 "scale.e"     :  0.465,
	#                 "p.y"         :  0.621,
	#                 "u.y"         : 68.084,
	#                 "sd.y"        :  8.626,
	#                 "shape.y"     :  0.057,
	#                 "p.hom"       :  0.853,
	#                 "u.hom"       : 11.101,
	#                 "sd.hom"      :  3.600,
	#                 "var.het"     : 10.916}

	sampleId = "apple_E12_L150_D80_K25"
	defaultParams = {"zp.copy.y"   :  3.000,
	                 "zp.copy.hom" :  3.000,
	                 "zp.copy.het" :  3.000,
	                 "p.e"         :  0.940,
	                 "shape.e"     :  3.000,
	                 "scale.e"     :  1.000,
	                 "p.y"         :  0.900,
	                 "u.y"         : 62.000,
	                 "sd.y"        : 16.309,
	                 "shape.y"     :  0.000,
	                 "p.hom"       :  0.800,
	                 "u.hom"       :  4.960,
	                 "sd.hom"      :  1.305,
	                 "var.het"     :  1.702}
	goodParams    = {"zp.copy.y"   :  2.047,
	                 "zp.copy.hom" :  3.390,
	                 "zp.copy.het" :  1.137,
	                 "p.e"         :  0.937,
	                 "shape.e"     :  0.114,
	                 "scale.e"     :  0.452,
	                 "p.y"         :  0.630,
	                 "u.y"         : 65.974,
	                 "sd.y"        :  8.666,
	                 "shape.y"     :  0.228,
	                 "p.hom"       :  0.818,
	                 "u.hom"       : 13.622,
	                 "sd.hom"      :  4.086,
	                 "var.het"     : 15.274}

	fitter = EnrichedHapDipFitter(path+"/"+sampleId+".mixed.kmer_dist")
	paramNames = fitter.paramNames

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
