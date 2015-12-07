#!/usr/bin/env bash
set -eu

# Create a sample of reads mixed from one chromosome and another and count the
# kmers in those reads.  Note that we discard the reads since all we want is
# the kmer distribution.

sampleId=$1
haploidChrom=$2
enhancement=$3   # ratio of reads_from_haploid / reads_from_diploid
readLen=$4
readDepth=$5
kmerSize=$6

# figure out weights for the sampler

haploidLen=`cat data/chrom.acgt_size \
  | awk '{ if ($1 == H) print $2; }' H=${haploidChrom}`
diploidLen=`cat data/test.diploid1.fa | fasta_lengths | awk '{ print $2; }'`
totalLen=$((haploidLen+diploidLen))
numReads=$((totalLen*readDepth/readLen))
haploidWeight=`echo ${enhancement} ${haploidLen} ${diploidLen} \
				 | awk '{ print $1/($2/$3) }'`
echo "haploidLen=${haploidLen} diploidLen=${diploidLen}"

# sample the reads

echo "sampling ${numReads} reads for ${sampleId}"
time simulate_reads --seed=${sampleId}.reads \
  ${haploidWeight}:data/hg19.${haploidChrom}.fa \
  .5:data/test.diploid1.fa \
  .5:data/test.diploid2.fa \
  ${numReads}x${readLen} \
  --prohibitN=5 \
  --noise=1% \
  --name=FAKE_[9]_{chrom}_{zstart}_{strand} \
  --width=none \
  --debug=input \
  --debug=choices \
  > farf.${sampleId}.mixed.reads.fa

# separate the two component sets so we can count kmers for them separately

echo "collecting haploid reads"
time cat farf.${sampleId}.mixed.reads.fa \
  | fasta_to_one_line --nosep \
  | awk '{ if ($1~H) { print $1;  print $2; }}' H=${haploidChrom} \
  > farf.${sampleId}.haploid_from_mixed.reads.fa

echo "collecting diploid reads"
time cat farf.${sampleId}.mixed.reads.fa \
  | fasta_to_one_line --nosep \
  | awk '{ if (!($1~H)) { print $1;  print $2; }}' H=${haploidChrom} \
  > farf.${sampleId}.diploid_from_mixed.reads.fa

# do the kmer counting ... the mixed set and the separated component sets

for ds in mixed haploid_from_mixed diploid_from_mixed ; do
	echo "counting kmers from ${sampleId}.${ds}"
	# abundance-min set very high, since all I need is the histogram
	time dsk -verbose 0 \
	  -file farf.${sampleId}.${ds}.reads.fa \
	  -abundance-min 200 \
	  -kmer-size ${kmerSize} \
	  -out farf.${sampleId}.${ds} \
	  -out-tmp temp/
	#
	h5dump --noindex \
		  --dataset=histogram/histogram \
		  farf.${sampleId}.${ds}.h5 \
	  | grep "^\ *[0-9]" \
	  | tr -d " " \
	  | tr -d "," \
	  | paste - - \
	  | awk '{ if ($1<500) print $0 }' \
	  > kmer_histograms/${sampleId}.${ds}.kmer_dist
	#
	rm farf.${sampleId}.${ds}.h5
	done

# cleanup temporary files

rm farf.${sampleId}.mixed.reads.fa 
rm farf.${sampleId}.haploid_from_mixed.reads.fa
rm farf.${sampleId}.diploid_from_mixed.reads.fa

