#! /usr/bin/env python2.7
"""
This script post-processes amplicon sequenced DBLa-tags basecalled using Multipass.
"""
import os, sys, time, copy
from optparse import OptionParser
import subprocess as Sp
from shlex import split as cmdsplit

# Constants
USEARCH = "/usr/local/bin/usearch_v5.2.32_i86osx32"
HMMSEARCH = "/usr/local/bin/hmmsearch_v3.1"
BLASTN = "/usr/local/bin/blastn_v2.2.25"
RESOURCES = "./resources"

# Script version
VERSION = "1.0.0"

#############################################################################################
#############################################################################################
# Command options

def build_parser():
	""" Builds command line parser """
	
	vers = "%%prog %s" %VERSION
	use = """usage: %prog [options] <basecalls_file>

<basecalls_file> File with '.basecalls'-extension produced by Multipass."""
	
	parser = OptionParser(usage=use, version=vers)
	
	parser.set_description("This script post-processes amplicon sequenced DBLa-tags basecalled using Multipass.")

	parser.add_option("-v", dest="v", action="store_true",
						help="print verbose information [False]",
						default=False)
	
	parser.add_option("-R", type="string", dest="R", metavar="RESULT_DIR",
						help="directory for result files [<basecall_file>.postprocess]",
						default="")

	parser.add_option("-i", dest="i", action="store_true",
						help="remove 3D7 sequences using BLASTN [False]",
						default=False)
	
	parser.add_option("-m", type="int", dest="m", metavar="METHOD",
						help="basecalling method. 0=Multipass, 1=Multipass_FRF  [1]",
						default=1)
	
	return parser

#############################################################################################
#############################################################################################
### Main action

def main():
	""" Main method """
	#########################
	global options
	### Parse commandline
	parser = build_parser()
	if len(sys.argv) < 2:
		parser.print_help(); sys.exit()
	(options, args) = parser.parse_args()
	if not os.path.isfile(args[0]):
		raise SystemExit("Error, could not find file: %s"% args[0])
	if not os.path.isdir(RESOURCES): 
		raise SystemExit("Error, could not find directory with resources: %s"% RESOURCES)
	fn = os.path.split(args[0])[1].strip(".basecalls")
	
	#########################
	
	### Create resultdir
	resultdir = options.R if options.R else args[0]+".postprocess"
	if not os.path.isdir(resultdir): 
		printv("Creating directory for results: %s"% resultdir)
		os.mkdir(resultdir)
	else: printv("Result directory: %s"% resultdir)
	resultPathFn = os.path.join(resultdir, fn)
	
	### Create fasta with most likely sequence for each flowgram cluster
	bcFastaFile = resultPathFn+'.bc.fas'
	with open(bcFastaFile, 'w') as fh:
		for (sid, seq) in getMLSeqs(args[0], fn, options.m):
			fh.write(">%s\n%s\n"%(sid, seq))
	
	### Identity clustering with seeds as output
	clusterSeedsFastaFile = resultPathFn+'.clu.fas'
	cmd = "%s --usersort --cluster %s --nofastalign --id 0.96 --seedsout %s --uc %s.clu.uc --sizein --sizeout"%(USEARCH, bcFastaFile, clusterSeedsFastaFile, resultPathFn)
	printv("Executing command: %s"%cmd); sys.stdout.flush()
	P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
	(sout, serr) = P.communicate()
	if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
	
	### Sort by size
	sortedSeedsFastaFile = resultPathFn+'.clus.fas'
	sortSize(clusterSeedsFastaFile, sortedSeedsFastaFile)
	
	### Chimera detection de-novo mode
	nonChimerasDenovo = resultPathFn+'.clus.nc.fas'
	cmd = "%s --uchime %s --uchimeout %s.clus.uchimeout --uchimealns %s.clus.uchimealns --chimeras %s.clus.ch.fas --nonchimeras %s"% \
		(USEARCH, sortedSeedsFastaFile, resultPathFn, resultPathFn, resultPathFn, nonChimerasDenovo)
	printv("Executing command: %s"%cmd); sys.stdout.flush()
	P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
	(sout, serr) = P.communicate()
	if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
	
	### Chimera detection database mode
	dbFile = resultPathFn+'.clus.nc.db'
	# Use clusters of at least size 2 as self (i.e. dont allow single sequence clusters to remove larger clusters)
	sortSize(nonChimerasDenovo, dbFile, mins=2)
	nonChimerasSelf = resultPathFn+'.clus.ncself.fas'
	# Remove chimeras by search against self  - increase sensitivity by setting --minh 0.2
	cmd = "%s --uchime %s --db %s --minh 0.2 --self --uchimeout %s.clus.uchimeoutself --uchimealns %s.clus.uchimealnsself --chimeras %s.clus.chself.fas --nonchimeras %s"% \
		(USEARCH, nonChimerasDenovo, dbFile, resultPathFn, resultPathFn, resultPathFn, nonChimerasSelf)
	printv("Executing command: %s"%cmd); sys.stdout.flush()
	P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
	(sout, serr) = P.communicate()
	if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
	
	### Remove least supported sequences
	trimmedFastaFile = resultPathFn+'.clust.fas'
	sortSize(nonChimerasSelf, trimmedFastaFile, minc=3)
	
	### Identify non-DBLa
	(sdic, slis) = fastaRead(trimmedFastaFile)
	nondbla = copy.deepcopy(slis)
	nonatagFastaFile = resultPathFn+'.nonatag.fas'
	cmd = "%s --cut_ga --domtblout %s.atag.hmmsearch -o /dev/null %s/atag.hmm %s"%(HMMSEARCH, resultPathFn, RESOURCES, trimmedFastaFile)
	printv("Executing command: %s"%cmd); sys.stdout.flush()
	P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
	(sout, serr) = P.communicate()
	if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
	with open(resultPathFn+".atag.hmmsearch") as fh:
		found = set([])
		for line in fh:
			if line.startswith("#"): continue
			ssl = line.strip().split()
			if ssl and (ssl[0] not in found):
				nondbla.remove(ssl[0])
				found.add(ssl[0])
	with open(nonatagFastaFile, 'w') as fh:
		for sid in nondbla:
			fh.write(">%s\n%s\n"%(sid, sdic[sid]))
	
	### Identify DBLb
	dblb = []
	btagFastaFile = resultPathFn+'.btag.fas'
	cmd = "%s --cut_ga --domtblout %s.btag.hmmsearch -o /dev/null %s/btag.hmm %s"%(HMMSEARCH, resultPathFn, RESOURCES, trimmedFastaFile)
	printv("Executing command: %s"%cmd); sys.stdout.flush()
	P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
	(sout, serr) = P.communicate()
	if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
	with open(resultPathFn+".btag.hmmsearch") as fh:
		found = set([])
		for line in fh:
			if line.startswith("#"): continue
			ssl = line.strip().split()
			if ssl and (ssl[0] not in found):
				dblb.append(ssl[0])
				found.add(ssl[0])
	with open(btagFastaFile, 'w') as fh:
		for sid in dblb:
			fh.write(">%s\n%s\n"%(sid, sdic[sid]))
	
	### Blast 3D7
	b3d7 = []
	if options.i:
		### Blast sequences against 3d7 DBLa-tags
		cmd = "%s -query %s -out /dev/stdout -db %s/3D7_dblatag/3D7dblatags.fas -evalue 1e-150 -perc_identity 96 -num_threads 2 -max_target_seqs 1 -outfmt '6 qseqid' "% \
			(BLASTN, trimmedFastaFile, RESOURCES)
		printv("Executing command: %s"%cmd); sys.stdout.flush()
		P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
		(sout,serr) = P.communicate()
		if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
		with open(resultPathFn+'.3d7atags.fas', 'w') as bfh:
			for line in sout.splitlines():
				sid = line.strip()
				if not sid: continue
				b3d7.append(sid)
				bfh.write(">%s\n%s\n"%(sid, sdic[sid]))
				
		
		### Blast sequences against 3d7 remaining genome
		cmd = "%s -strand 'plus' -query %s -out /dev/stdout -db %s/3D7_nondblatag/PlasmoDB-9.3_Pfalciparum3D7_Genome_mar2013.allfwd.fas -evalue 1e-50 -perc_identity 96 -num_threads 2 -max_target_seqs 1 -outfmt '6 qseqid' "% \
			(BLASTN, trimmedFastaFile, RESOURCES)
		printv("Executing command: %s"%cmd); sys.stdout.flush()
		P = Sp.Popen(cmdsplit(cmd), stdout=Sp.PIPE, stdin=Sp.PIPE, stderr=Sp.PIPE)
		(sout,serr) = P.communicate()
		if P.wait() != 0: raise SystemExit("Problem executing: %s\n%s"%(cmd, serr)) 
		with open(resultPathFn+'.3d7nonatags.fas', 'w') as bfh:
			for line in sout.splitlines():
				sid = line.strip()
				if not sid: continue
				b3d7.append(sid)
				bfh.write(">%s\n%s\n"%(sid, sdic[sid]))
	
	### Write clean DBLa-tags
	atagFastaFile = resultPathFn+'.clean.fas'
	with open(atagFastaFile, 'w') as fh:
		for sid in slis:
			if sid not in (dblb+nondbla+b3d7):
				fh.write(">%s\n%s\n"%(sid, sdic[sid]))
	

#############################################################################################
#############################################################################################
### Sort fasta file by cluster size

def sortSize(inFile, outFile, minc=1, mins=1):
	""" Sort fasta file by cluster size """
	(sdic, slis) = fastaRead(inFile)
	lol = []
	for sid in slis:
		remain = sid
		(v, remain) = remain.rsplit("_c",1)
		(c, s) = remain.rsplit(";size=")
		if (int(c) >= minc) and (int(s) >= mins):
			lol.append((int(s), int(c), sid))
	lol.sort(reverse = True)
	with open(outFile, 'w') as outFileHand:
		for (s, c, sid) in lol:
			outFileHand.write(">%s\n%s\n"%(sid, sdic[sid]))

#############################################################################################
#############################################################################################
### Read fasta file

def fastaRead(fasFile):
	""" Read fasta file """
	seqDict = {}
	seqList = []
	with open(fasFile) as fasFileHand:
		for line in fasFileHand:
			sline=line.strip()
			if sline=='' or sline[0]=='#': continue
			elif sline[0]=='>':
				sid=sline[1:]
				seqDict[sid]=""
				seqList.append(sid)
			else:
				seqDict[sid] += sline.upper()
	return(seqDict, seqList)

#############################################################################################
#############################################################################################
### Get most likely sequence from each flowgram cluster

def getMLSeqs(mpFileName, fn, meth):
	""" Get most likely sequence from each flowgram cluster """
	mLSeqs = []
	for line in open(mpFileName):
		if line.startswith(">Cluster"): 
			if mLSeqs:
				if len(mLSeqs) < 10: raise SystemExit("Sorry, this method requires at least 10 most likely sequences per cluster.")
				mLSeq = max(mLSeqs[:10])
				sid = "%s_%s_c%s;size=%s"%(fn, clusterInfo[2].replace("_","."), cSize, cSize)
				yield(sid, mLSeq[1])
				mLSeqs = []
			clusterInfo = line.strip().split()
			cSize = int(clusterInfo[3])
			readML = False
		elif line.startswith(">>>"): 
			readML = True
		elif readML:
			ssline = line.strip().split()
			if len(ssline) == 3:
				if (meth == 0) or (cSize == 1): mLSeqs.append([float(ssline[1]), ssline[2]])
				elif (meth == 1): mLSeqs.append([float(ssline[1]) * Pcbs_frf(ssline[2]), ssline[2]])
				else: raise SystemExit("Illegal method value: %s"%(meth))
	if mLSeqs:
		mLSeq = max(mLSeqs)
		sid = "%s_%s_c%s;size=%s"%(fn, clusterInfo[2].replace("_","."), cSize, cSize)
		yield(sid, mLSeq[1])

#############################################################################################
#############################################################################################
### Probability that the sequence is correctly basecalled given the presence/absence
### of a full length reading frame

def Pcbs_frf(seq):
	if findORF(seq): return(0.489)
	else: return(2.58e-4)

#############################################################################################
#############################################################################################
### Find ORF in forward reading frame 1, 2, and 3

stdDnaCode = {'---':'-','TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
stopCodon = ("TAG", "TGA", "TAA")

def findORF(s):
	for k in range(3):
		orf = []
		for i in range(k, len(s), 3):
			cod = s[i:i+3].upper()
			if (len(cod) != 3): continue
			if cod in stopCodon: break
			else: orf.append(stdDnaCode.get(cod,'X'))
		else:
			return("".join(orf))

#############################################################################################
#############################################################################################
### Print verbose

def printv(txt):
	if options.v: print(txt)

#############################################################################################
#############################################################################################
### Main

if __name__ == "__main__":
    main()



