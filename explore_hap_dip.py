#!/usr/bin/env python
"""
Explore setting of initial parameters for hap-dip model, by seeing how far any
of the parameters can vary from the "right" answer before we lose convergence.
"""

from sys  import argv,stderr
from math import sqrt
from kmervature import EnrichedHapDipFitter,params_to_text


def main():
	assert (len(argv) == 1), "give me no arguments"
	path = "kmer_histograms"

	sampleId = "mixedB"
	fitter = EnrichedHapDipFitter(path+"/"+sampleId+".mixed.kmer_dist")
	paramNames = fitter.paramNames

	defaultParams = {"zp.copy.y"   :  3.000,
	                 "zp.copy.hom" :  3.000,
	                 "zp.copy.het" :  3.000,
	                 "p.e"         :  0.942,
	                 "shape.e"     :  3.000,
	                 "scale.e"     :  1.000,
	                 "p.y"         :  0.900,
	                 "u.y"         : 64.000,
	                 "sd.y"        : 14.826,
	                 "shape.y"     :  0.000,
	                 "p.hom"       :  0.800,
	                 "u.hom"       :  5.120,
	                 "sd.hom"      :  1.186,
	                 "var.het"     :  1.407}

	goodParams    = {"zp.copy.y"   :  2.042,
	                 "zp.copy.hom" :  3.157,
	                 "zp.copy.het" : 17.795,
	                 "p.e"         :  0.935,
	                 "shape.e"     :  0.096,
	                 "scale.e"     :  0.465,
	                 "p.y"         :  0.621,
	                 "u.y"         : 68.084,
	                 "sd.y"        :  8.626,
	                 "shape.y"     :  0.057,
	                 "p.hom"       :  0.853,
	                 "u.hom"       : 11.101,
	                 "sd.hom"      :  3.600,
	                 "var.het"     : 10.916}

	numSteps = 10
	for (paramIx,name) in enumerate(paramNames):
		if (paramIx != 0): print
		for step in xrange(1,numSteps+1):
			print "=== param %d of %s (\"%s\") step %d of %s ===" \
			    % (1+paramIx,len(paramNames),name,step,numSteps)

			initParams = dict(goodParams)
			initParams[name] += step*(defaultParams[name]-goodParams[name])/numSteps
			fitter.set_params(initParams)
			fitParams = fitter.fit()
			if (fitParams == None):
				print params_to_text(paramNames,initParams,prefix="init:")
				print "(failure or non-convergence)"
				print "... return code ..."
				print fitter.retCode
				print "... stdout ..."
				print fitter.stdout
				print "... stderr ..."
				print fitter.stderr
				continue

			print params_to_text(paramNames,initParams,fitParams,prefix="init:",prefix2="cvrg:")
			fitParams = params_to_float(fitParams)
			distance = vector_distance(fitParams,goodParams)
			print "dGood: %.8f" % distance


def params_to_float(params):
	return {name:float(params[name]) for name in params}


def vector_distance(vector1,vector2):
	if ([name for name in vector1 if (name not in vector2)] != []): raise valueError
	if ([name for name in vector2 if (name not in vector1)] != []): raise valueError
	return sqrt(sum([((vector1[name]-vector2[name])**2) for name in vector1]))


if __name__ == "__main__": main()
