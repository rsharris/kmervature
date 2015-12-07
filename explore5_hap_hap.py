#!/usr/bin/env python
"""
Without cheating, use estimates from a haploid fit to the mixed data to
contribute to the initial vector for the hap-hap fit.
"""

from sys import argv,stderr
from kmervature import HaploidFitter,DiploidFitter,EnrichedHapHapFitter, \
                       params_to_text,params_from_text


def main():
	assert (len(argv) == 2), "need the sampleID and nothing else"
	sampleId = argv[1]
	explainFailure = True
	path = "kmer_histograms"

	print sampleId

	# perform haploid fit to the sample (ignoring thge diploid component)

	hFitter = HaploidFitter(path+"/"+sampleId+".mixed.kmer_dist")
	hParamNames = hFitter.paramNames

	hFitParams = hFitter.fit()
	if (hFitParams == None):
		print >>stderr, "haploid: failure or non-convergence"
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

	# ask for default values for the enrichment model

	hhFitter = EnrichedHapHapFitter(path+"/"+sampleId+".mixed.kmer_dist")
	hhParamNames = hhFitter.paramNames

	hhDefaultParams = hhFitter.default_params()
	if (hhDefaultParams == None):
		print >>stderr, "hap-hap: failed to get default params"
		print "(hap-hap: failed to get default params)"
		if (explainFailure):
			print "... return code ..."
			print hhFitter.retCode
			print "... stdout ..."
			print hhFitter.stdout
			print "... stderr ..."
			print hhFitter.stderr
	else:
		print params_to_text(hhParamNames,hhDefaultParams,prefix="dflt.haphap:")

	assert (hFitParams != None) and (hhDefaultParams != None), \
	       "(no point in trying to fit the hap-hap model)"

	# create an initial vector for the enrichment model, borrowing some
	# elements from the haploid model fit

	hhInitParams = dict(hhDefaultParams)
	hhInitParams["zp.copy.y"] = hFitParams["zp.copy"]
	hhInitParams["p.e"      ] = hFitParams["p.e"]
	hhInitParams["shape.e"  ] = hFitParams["shape.e"]
	hhInitParams["scale.e"  ] = hFitParams["scale.e"]
	hhInitParams["u.y"      ] = hFitParams["u.v"]
	hhInitParams["sd.y"     ] = hFitParams["sd.v"]
	hhInitParams["shape.y"  ] = hFitParams["shape.v"]

	pAuto = 1 - float(hhInitParams["p.y"])
	hhInitParams["u.auto"   ] =         pAuto * float(hhInitParams["u.y"])
	hhInitParams["sd.auto"  ] = sdHom = pAuto * float(hhInitParams["sd.y"])

	# perform hap-hap fit to the mixed components

	hhFitParams = hhFitter.fit(hhInitParams)
	if (hhFitParams == None):
		print >>stderr, "hap-hap: failure or non-convergence"
		print "(hap-hap: failure or non-convergence)"
		print params_to_text(hhParamNames,hhInitParams,prefix="smart.haphap:")
		if (explainFailure):
			print "... return code ..."
			print hhFitter.retCode
			print "... stdout ..."
			print hhFitter.stdout
			print "... stderr ..."
			print hhFitter.stderr
	else:
		print params_to_text(hhParamNames,hhInitParams,hhFitParams,
		                     prefix="smart.haphap:",prefix2="cvrg.haphap:")


if __name__ == "__main__": main()
