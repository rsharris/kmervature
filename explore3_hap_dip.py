#!/usr/bin/env python
"""
Explore setting of initial parameters for hap-dip model, by trying fits for
both the haploid and diploid components (these are known *only* for simulated
datasets), and combining results from those fits into an inital vector for the
hap-dip model.

This is "cheating", since for real datasets we will not know the components.
But this provides a working answer for simulated datasets which can be used in
experimenting with variation of parameters.
"""

from sys    import argv,stderr
from math   import sqrt
from random import seed as random_seed,random as unit_random
from kmervature import HaploidFitter,DiploidFitter,EnrichedHapDipFitter,params_to_text


def main():
	assert (len(argv) == 2), "need the sampleID and nothing else"
	sampleId = argv[1]
	saveConvergence = True
	explainFailure = True
	path = "kmer_histograms"

	# clear the convergence file, in case we have a failure (we don't want
	# previous results to leak through)

	if (saveConvergence):
		f = file(path+"/"+sampleId+".mixed.fit","wt")
		f.close()

	# perform haploid fit to the haploid component

	hFitter = HaploidFitter(path+"/"+sampleId+".haploid_from_mixed.kmer_dist")
	hParamNames = hFitter.paramNames

	hFitParams = hFitter.fit()
	if (hFitParams == None):
		print "(haploid: failure or non-convergence)"
		if (explainFailure):
			print "... return code ..."
			print hFitter.retCode
			print "... stdout ..."
			print hFitter.stdout
			print "... stderr ..."
			print hFitter.stderr
	else:
		print params_to_text(hParamNames,hFitParams,prefix="cvrg.haploid:")

	# perform diploid fit to the diploid component

	dFitter = DiploidFitter(path+"/"+sampleId+".diploid_from_mixed.kmer_dist")
	dParamNames = dFitter.paramNames

	dFitParams = dFitter.fit()
	if (dFitParams == None):
		print "(diploid: failure or non-convergence)"
		if (explainFailure):
			print "... return code ..."
			print dFitter.retCode
			print "... stdout ..."
			print dFitter.stdout
			print "... stderr ..."
			print dFitter.stderr
	else:
		print params_to_text(dParamNames,dFitParams,prefix="cvrg.diploid:")

	# create an initial vector for the enrichment model, combining elements
	# from the component fits with the usual defaults

	hdFitter = EnrichedHapDipFitter(path+"/"+sampleId+".mixed.kmer_dist")
	hdParamNames = hdFitter.paramNames

	hdDefaultParams = hdFitter.default_params()
	if (hdDefaultParams == None):
		print "(hap-dip: failed to get default params)"
		if (explainFailure):
			print "... return code ..."
			print hdFitter.retCode
			print "... stdout ..."
			print hdFitter.stdout
			print "... stderr ..."
			print hdFitter.stderr
	else:
		print params_to_text(hdParamNames,hdDefaultParams,prefix="dflt.hapdip:")

	assert (hFitParams != None) and (dFitParams != None) and (hdDefaultParams != None), \
	       "(no point in trying to fit the hap-dip model)"

	hdInitParams = {}
	hdInitParams["zp.copy.y"  ] = hFitParams["zp.copy"]
	hdInitParams["zp.copy.hom"] = dFitParams["zp.copy"]
	hdInitParams["zp.copy.het"] = dFitParams["zp.copy.het"]
	hdInitParams["p.e"        ] = hFitParams["p.e"]
	hdInitParams["shape.e"    ] = hFitParams["shape.e"]
	hdInitParams["scale.e"    ] = hFitParams["scale.e"]
	hdInitParams["p.y"        ] = hdDefaultParams["p.y"]
	hdInitParams["u.y"        ] = hFitParams["u.v"]
	hdInitParams["sd.y"       ] = hFitParams["sd.v"]
	hdInitParams["shape.y"    ] = hFitParams["shape.v"]
	hdInitParams["p.hom"      ] = 1 - float(dFitParams["p.d"])
	hdInitParams["u.hom"      ] = dFitParams["u.v"]
	hdInitParams["sd.hom"     ] = dFitParams["sd.v"]
	hdInitParams["var.het"    ] = dFitParams["var.w"]

	# perform hap-dip fit to the mixed components

	hdFitParams = hdFitter.fit(hdInitParams)
	if (hdFitParams == None):
		print "(hap-dip: failure or non-convergence)"
		print params_to_text(hdParamNames,hdInitParams,prefix="init.hapdip:")
		if (explainFailure):
			print "... return code ..."
			print hdFitter.retCode
			print "... stdout ..."
			print hdFitter.stdout
			print "... stderr ..."
			print hdFitter.stderr
	else:
		print params_to_text(hdParamNames,hdInitParams,hdFitParams,
		                     prefix="init.hapdip:",prefix2="cvrg.hapdip:")

	# write the convergence file

	if (saveConvergence):
		f = file(path+"/"+sampleId+".mixed.fit","wt")
		print >>f, params_to_text(hdParamNames,hdFitParams)
		f.close()


if __name__ == "__main__": main()
