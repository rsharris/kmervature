#!/usr/bin/env python
"""
Without cheating, use values from a haploid fit to the mixed data to contribute
to the initial vector for the hap-dip fit.  If that fails, try to determine how
far away from success it is.
"""

from sys import argv,stderr
from kmervature import HaploidFitter,DiploidFitter,EnrichedHapDipFitter, \
                       params_to_text,params_from_text


def main():
	assert (len(argv) == 2), "need the sampleID and nothing else"
	sampleId = argv[1]
	explainFailure = True
	path = "kmer_histograms"
	paramsToCheat = {}
	#paramsToCheat["p.y"] = 0.5   # 1.0 means "cheat completely; 0 means not at all"
	paramsToFudge = {}
	#paramsToFudge["shape.e"] = 0.55
	#paramsToFudge["scale.e"] = 1.40

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

	hdFitter = EnrichedHapDipFitter(path+"/"+sampleId+".mixed.kmer_dist")
	hdParamNames = hdFitter.paramNames

	hdDefaultParams = hdFitter.default_params()
	if (hdDefaultParams == None):
		print >>stderr, "hap-dip: failed to get default params"
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

	assert (hFitParams != None) and (hdDefaultParams != None), \
	       "(no point in trying to fit the hap-dip model)"

	# read the sample's "cheat" parameters for comparison (usually produced by
	# explore3_hap_dip)

	fitFilename = path+"/"+sampleId+".mixed.fit"

	f = file(fitFilename,"rt")
	hdCheatParams = params_from_text([line for line in f])
	f.close()

	for name in hdDefaultParams:
		assert (name in hdCheatParams), \
		       "parameter \"%s\" missing from %s" % (name,fitFilename)

	for name in hdCheatParams:
		assert (name in hdDefaultParams), \
		       "extra parameter \"%s\" in %s" % (name,fitFilename)

	# create an initial vector for the enrichment model, borrowing some
	# elements from the haploid model fit

	hdInitParams = dict(hdDefaultParams)
	hdInitParams["zp.copy.y"] = hFitParams["zp.copy"]
	hdInitParams["p.e"      ] = hFitParams["p.e"]
	hdInitParams["shape.e"  ] = hFitParams["shape.e"]
	hdInitParams["scale.e"  ] = hFitParams["scale.e"]
	hdInitParams["u.y"      ] = hFitParams["u.v"]
	hdInitParams["sd.y"     ] = hFitParams["sd.v"]
	hdInitParams["shape.y"  ] = hFitParams["shape.v"]

	for name in paramsToCheat:
		param = float(hdInitParams[name])
		param += paramsToCheat[name] * (float(hdCheatParams[name]) - param)
		hdInitParams[name] = param
	for name in paramsToFudge:
		hdInitParams[name] = float(hFitParams[name]) * paramsToFudge[name]

	pAuto = 1 - float(hdInitParams["p.y"])
	pHom  =     float(hdInitParams["p.hom"])
	hdInitParams["u.hom"    ] =         pAuto * pHom * float(hdInitParams["u.y"])
	hdInitParams["sd.hom"   ] = sdHom = pAuto * pHom * float(hdInitParams["sd.y"])
	hdInitParams["var.het"  ] = sdHom * sdHom

	# perform hap-dip fit to the mixed components

	hdFitParams = hdFitter.fit(hdInitParams)
	if (hdFitParams == None):
		print >>stderr, "hap-dip: failure or non-convergence"
		print "(hap-dip: failure or non-convergence)"
		print params_to_text(hdParamNames,hdInitParams,hdCheatParams,
		                     prefix="smart.hapdip:",prefix2="cheat.hapdip:")
		if (explainFailure):
			print "... return code ..."
			print hdFitter.retCode
			print "... stdout ..."
			print hdFitter.stdout
			print "... stderr ..."
			print hdFitter.stderr
	else:
		print params_to_text(hdParamNames,hdInitParams,hdFitParams,
		                     prefix="smart.hapdip:",prefix2="cvrg.hapdip:")
		print params_to_text(hdParamNames,hdCheatParams,prefix="cheat.hapdip:")

	# if convergence failed, try moving the initial parameters toward the
	# cheat parameters in small steps until we get convergence
	# $$$ a binary search would be "better"

	numSteps = 100
	step = 0

	while (hdFitParams == None):
		step += 1
		if (step == numSteps): break
		print >>stderr, "step %d" % step

		hdStepParams = {}
		for name in hdInitParams:
			if (name in ["u.hom","sd.hom","var.het"]): continue
			param = float(hdInitParams[name])
			param += (step * (float(hdCheatParams[name]) - param)) / numSteps
			hdStepParams[name] = param

		pAuto = 1 - float(hdStepParams["p.y"])
		pHom  =     float(hdStepParams["p.hom"])
		hdStepParams["u.hom"    ] =         pAuto * pHom * float(hdStepParams["u.y"])
		hdStepParams["sd.hom"   ] = sdHom = pAuto * pHom * float(hdStepParams["sd.y"])
		hdStepParams["var.het"  ] = sdHom * sdHom

		hdFitParams = hdFitter.fit(hdStepParams)
		if (hdFitParams == None):
			print params_to_text(hdParamNames,hdStepParams,
			                     prefix="step[%d].hapdip:" % step)
			#if (explainFailure):
			#	print "... return code ..."
			#	print hdFitter.retCode
			#	print "... stdout ..."
			#	print hdFitter.stdout
			#	print "... stderr ..."
			#	print hdFitter.stderr
		else:
			print params_to_text(hdParamNames,hdStepParams,hdFitParams,
			                     prefix="step[%d].hapdip:" % step,prefix2="cvrg.hapdip:")


if __name__ == "__main__": main()
