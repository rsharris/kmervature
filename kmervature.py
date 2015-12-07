#!/usr/bin/env python
"""
kmervature-- model fitting to kmer abundance histograms
Copyright (C) 2015 Bob Harris

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

===

Wrappers for curve-fitting functions in kmervature.r.

A simple test program is included, as main().
"""

from sys import argv,stderr
import subprocess


# HaploidFitter class

class HaploidFitter(object):

	def __init__(self,histFile,initParams=None,rscriptPath=None):
		self.histFile   = histFile
		self.modelType  = "haploid"
		self.paramNames = ["zp.copy",
		                   "p.e","shape.e","scale.e",
		                   "u.v","sd.v","shape.v"]
		self.common_init(initParams=initParams,rscriptPath=rscriptPath)

	def common_init(self,initParams=None,rscriptPath=None):
		"""
		not expected to be called externally
		"""
		self.fpPrecision = 8

		if   (rscriptPath == None):       self.rscriptPath = ""
		elif (rscriptPath.endswith("/")): self.rscriptPath = rscriptPath
		else:                             self.rscriptPath = rscriptPath + "/"

		self.initParams = {}
		for name in self.paramNames:
			self.initParams[name] = None
		if (initParams != None):
			self.set_params(initParams)

		self.fitParams = {}
		self.retCode = self.stdout = self.stderr = None
		self.debug = []

	def set_params(self,params):
		"""
		changes some or all of the initial parameters to be used by fit();
		params can be named values in a dict, or they can be a list with values
		in the same order as in paramNames
		"""
		if (type(params) == dict):
			for name in self.paramNames:
				if (name not in params): continue
				self.initParams[name] = params[name]
		elif (type(params) == list):
			if (len(params) != len(self.paramNames)): raise ValueErr
			for (ix,name) in enumerate(self.paramNames):
				self.initParams[name] = params[ix]
		else:
			raise ValueErr

	def default_params(self):
		"""
		returns a dict mapping parameter names to values (as strings);  None
		if failure
		"""

		# create a short one-line R program;  note that here each command is
		# a separate string in a list, and separating semi-colons will be added
		# when we go to run it

		rCommand =  []
		rCommand += ["source('%skmervature.r')" % self.rscriptPath]
		rCommand += ["hist.file='%s'" % self.histFile]
		rCommand += ["cov = read.coverage(hist.file)"]
		rCommand += ["fit = kmer.histogram.fit('%s',cov,max.copy=3,performFit=F)" % self.modelType]
		rCommand += ["cat('params.init=')"]
		rCommand += ["for (i in 1:length(fit$params.init)) { cat(sprintf(' %%.%df',fit$params.init[i])) }" % self.fpPrecision]
		rCommand += ["cat('\r')"]

		# run that R program

		(retCode,out,err) = self.rscript(["-e","; ".join(rCommand)])
		self.retCode = retCode
		self.stdout  = out
		self.stderr  = err

		if ("show output" in self.debug):
			print >>stderr, "=== return code ==="
			print >>stderr, retCode
			print >>stderr, "=== stdout ==="
			print >>stderr, out
			print >>stderr, "=== stderr ==="
			print >>stderr, err

		if (retCode != 0): return None

		return self.parse_for_params(out,varName="params.init")

	def fit(self,params=None):
		"""
		returns a dict mapping parameter names to values (as strings);  None
		if failure
		"""

		if (params != None):
			self.set_params(params)

		# convert initial params, if they exist, to strings;  note that while
		# we allow the caller to use set_params to set different parameters
		# willy-nilly, here we require either that *all* parameters have been
		# given a value, or *none*;  values that were provided as strings are
		# not altered here;  and empty string and None indicate parameters
		# that have not been set

		initParams = []
		setParam = None
		unsetParam = None
		for name in self.paramNames:
			paramVal = self.initParams[name]
			if (paramVal == None):
				unsetParam = name
				continue
			elif (type(paramVal) == str):
				if (paramVal == ""):
					unsetParam = name
					continue
			else:
				paramVal = "%.*f" % (self.fpPrecision,paramVal)
			initParams += [paramVal]
			setParam = name

		if (setParam == None):
			initParams = None
		elif (unsetParam != None):
			raise RuntimeError("for %s, %s has been given an initial value but %s has not" \
			                 % (self.modelType,setParam,unsetParam))

		# create a short one-line R program;  note that here each command is
		# a separate string in a list, and separating semi-colons will be added
		# when we go to run it

		rCommand =  []
		rCommand += ["source('%skmervature.r')" % self.rscriptPath]
		rCommand += ["hist.file='%s'" % self.histFile]
		rCommand += ["cov = read.coverage(hist.file)"]
		if (initParams != None):
			rCommand += ["init.params=c(%s)" % ",".join(initParams)]
			rCommand += ["fit = kmer.histogram.fit('%s',cov,max.copy=3,params.init=init.params)" % self.modelType]
		else:
			rCommand += ["fit = kmer.histogram.fit('%s',cov,max.copy=3)" % self.modelType]
		rCommand += ["cat('fit.params=')"]
		rCommand += ["for (i in 1:length(fit$params)) { cat(sprintf(' %%.%df',fit$params[i])) }" % self.fpPrecision]
		rCommand += ["cat('\r')"]

		# run that R program;  note that we neuter the fitParams result, so
		# that a result from an earlier run won't leak through if this run
		# fails

		self.fitParams = {}
		(retCode,out,err) = self.rscript(["-e","; ".join(rCommand)])
		self.retCode = retCode
		self.stdout  = out
		self.stderr  = err

		if ("show output" in self.debug):
			print >>stderr, "=== return code ==="
			print >>stderr, retCode
			print >>stderr, "=== stdout ==="
			print >>stderr, out
			print >>stderr, "=== stderr ==="
			print >>stderr, err

		if (retCode != 0): return None

		self.fitParams = self.parse_for_params(out)
		return self.fitParams

	def parse_for_params(self,txt,varName="fit.params="):
		"""
		not expected to be called externally
		"""
		if (not varName.endswith("=")): varName += "="

		paramTxt = None
		for line in txt.split("\n"):
			if (line.startswith(varName)):
				paramTxt = line
				break
		if (paramTxt == None):
			return None

		paramTxt = paramTxt.split()[1:]
		if (len(paramTxt) != len(self.paramNames)):
			return None

		fitParams = {}
		for (ix,name) in enumerate(self.paramNames):
			fitParams[name] = paramTxt[ix]

		return fitParams

	def rscript(self,rCommand):
		"""
		not expected to be called externally
		"""
		return self.run(["Rscript","--no-init-file"] + rCommand)

	def run(self,command):
		"""
		not expected to be called externally
		"""
		process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		(out,err) = process.communicate()
		retCode = process.wait()
		return (retCode,out,err)


# other Fitter classes

class DiploidFitter(HaploidFitter):

	def __init__(self,histFile,initParams=None,rscriptPath=None):
		self.histFile    = histFile
		self.modelType   = "diploid"
		self.paramNames = ["zp.copy",
		                   "p.e","shape.e","scale.e",
		                   "u.v","sd.v",
		                   "p.d","var.w","zp.copy.het"]
		self.common_init(initParams=initParams,rscriptPath=rscriptPath)


class EnrichedHapHapFitter(HaploidFitter):

	def __init__(self,histFile,initParams=None,rscriptPath=None):
		self.histFile    = histFile
		self.modelType   = "enriched-hap-hap"
		self.paramNames = ["zp.copy.y","zp.copy.auto",
		                   "p.e","shape.e","scale.e",
		                   "p.y","u.y","sd.y","shape.y",
		                   "u.auto","sd.auto"]
		self.common_init(initParams=initParams,rscriptPath=rscriptPath)


class EnrichedHapDipFitter(HaploidFitter):

	def __init__(self,histFile,initParams=None,rscriptPath=None):
		self.histFile   = histFile
		self.modelType  = "enriched-hap-dip"
		self.paramNames = ["zp.copy.y","zp.copy.hom","zp.copy.het",
		                   "p.e","shape.e","scale.e",
		                   "p.y","u.y","sd.y","shape.y",
		                   "p.hom","u.hom","sd.hom","var.het"]
		self.common_init(initParams=initParams,rscriptPath=rscriptPath)


# some classless support functions

def params_to_text(hdParamNames,params,params2=None,prefix=None,prefix2=None):
	fieldW = {}
	for name in hdParamNames:
		fieldW[name] = max(len(name),len(str(params[name])))
		if (params2 != None):
			fieldW[name] = max(fieldW[name],len(str(params2[name])))

	prefixW = 0
	if (prefix  != None): prefixW = len(prefix)+1
	if (prefix2 != None): prefixW = max(prefixW,len(prefix2)+1)

	lines = []
	if (prefixW == 0):
		lines += [" ".join(["%-*s" % (fieldW[name],name)          for name in hdParamNames])]
		lines += [" ".join(["%-*s" % (fieldW[name],params[name])  for name in hdParamNames])]
		if (params2 != None):
			lines += [" ".join(["%-*s" % (fieldW[name],params2[name]) for name in hdParamNames])]
	else:
		hdrPrefix = "." + (" " * (prefixW-1))
		if (prefix  == None): prefix  = " " * prefixW
		else:                 prefix  = "%-*s" % (prefixW-1,prefix)
		if (prefix2 == None): prefix2 = " " * prefixW
		else:                 prefix2 = "%-*s" % (prefixW-1,prefix2)

		lines += [hdrPrefix     + " ".join(["%-*s" % (fieldW[name],name)          for name in hdParamNames])]
		lines += [prefix  + " " + " ".join(["%-*s" % (fieldW[name],params[name])  for name in hdParamNames])]
		if (params2 != None):
			lines += [prefix2 + " " + " ".join(["%-*s" % (fieldW[name],params2[name]) for name in hdParamNames])]

	return "\n".join(lines)


def params_from_text(txt):
	if (type(txt) == list): lines = list(txt)
	else:                   lines = txt.split("\n")
	if (len(lines) != 2): raise ValueError("%d lines of text (expected exactly 2)" % len(lines))

	names = lines[0].split()
	if (len(names) < 1): raise ValueError("first line contains no parameter names")

	vals = lines[1].split()
	if (len(vals) != len(names)): raise ValueError("%d parameter names but %d values" % (len(names),len(vals)))

	nameToVal = {}
	for (ix,name) in enumerate(names):
		if (name in nameToVal):
			raise ValueError("parameter name \"%s\" appears more than once" % name)
		nameToVal[name] = vals[ix]

	return nameToVal


# simple test program

def main():
	assert (len(argv) == 2), "need model type and nothing else"
	modelType = argv[1]

	path = "kmer_histograms"

	if (modelType == "haploid"):
		fitter = HaploidFitter("%s/mixedB.haploid_from_mixed.kmer_dist" % path)

	elif (modelType == "diploid"):
		fitter = DiploidFitter("%s/mixedB.diploid_from_mixed.kmer_dist" % path)

	elif (modelType == "enriched-hap-hap"):
		fitter = EnrichedHapHapFitter("%s/mixedB.mixed.kmer_dist" % path)
		fitter.set_params({"zp.copy.y"    :  2.001,
		                   "zp.copy.auto" :  3.577,
		                   "p.e"          :  0.956,
		                   "shape.e"      :  0.209,
		                   "scale.e"      :  0.431,
		                   "p.y"          :  0.900,
		                   "u.y"          : 70.386,
		                   "sd.y"         :  8.986,
		                   "shape.y"      : -0.285,
		                   "u.auto"       : 11.148,
		                   "sd.auto"      :  3.602})

	elif (modelType == "enriched-hap-dip"):
		fitter = EnrichedHapDipFitter("%s/mixedB.mixed.kmer_dist" % path)
		fitter.set_params({"zp.copy.y"   :  2.042,
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
		                   "var.het"     : 10.916})

	else:
		assert (False), "unkown model type: \"%s\"" % modelType

	#fitter.debug += ["show output"]

	print "=== defaults ==="
	defaultParams = fitter.default_params()
	if (defaultParams == None):
		print "defaults failed!"
	else:
		nameW = max(len(name) for name in fitter.paramNames)
		for name in fitter.paramNames:
			print "%-*s = %s" % (nameW,name,defaultParams[name])

	print "=== convergence ==="
	fitParams = fitter.fit()
	if (fitParams == None):
		print "did not converge"
	else:
		nameW = max(len(name) for name in fitter.paramNames)
		for name in fitter.paramNames:
			print "%-*s = %s" % (nameW,name,fitParams[name])


if __name__ == "__main__": main()
