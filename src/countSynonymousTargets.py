#!/usr/bin/python

## countSynonymousTargets.py by Rohan Maddamsetti.
## This script takes a genbank file (i.e. REL606.gbk),
## and produces a csv file that counts the number of possible sites where a
## synonymous G->A, G->T, etc. mutation could occur within a protein-coding
## locus.

import sys
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from csv import writer

def mutate_codon(codon, before, after, trans_table):
	"""Make a list of codons with single substititions of letters from->to.
	for example: mutate_codon('GTG','G','A', trsl_table) 
	produces ['ATG', 'GTA'], then count how many of these mutations 
	are synonymous."""
	amino_acid = codon.translate(table=trans_table)
	mutated = []
	for i,letter in enumerate(codon):
		if letter == before:
			codon_list = list(codon)
			codon_list[i] = after
			mutated.append("".join(codon_list))
	translated = [SeqIO.Seq(x).translate(table=trans_table) for x in mutated]
	synonymous_changes = sum([1 for y in translated if str(y) == str(amino_acid)])
	return synonymous_changes

def codons(s): 
	"""Return list of codons for the dna string.""" 
	end = len(s) - (len(s) % 3) - 1 
	codons = [s[i:i+3] for i in range(0, end, 3)] 
	return codons 

def main():
	input_gbk = None
	if (len(sys.argv) == 2):
		input_gbk = sys.argv[1]
	else:
		sys.exit("ERROR: Wrong number of arguments provided!")
	
	## Save the data in a dict of dicts.
	## Key: locus_tag, Value: the following dict of data;
	## Key: 'possible.G.to.A' (etc.), Value: the number of sites.
	locus_dict = {}

	## Iterate through all CDS in the genbank file.
	genome = SeqIO.parse(input_gbk, "genbank").next()
	for feature in genome.features:
		if (feature.type != "CDS"):
			continue
		## Get the locus_tag.
		locus_tag = feature.qualifiers["locus_tag"][0]
		## Get the name of the gene.
		try:
			gene_name = feature.qualifiers["gene"][0]
		except KeyError:
			gene_name = locus_tag
		dna_sequence = feature.extract(genome.seq)
		trans_table = feature.qualifiers["transl_table"][0]
		## check that the translation is correct;
		## the exceptions are pseudogenes, or contain selenocysteine.
		try: 
			protein = dna_sequence.translate(table=trans_table,cds=True)
		except TranslationError:
			continue
		## Initialize the dict of data.
		sitecounts = dict.fromkeys(["possible.G.to.A", "possible.G.to.T", 
									"possible.G.to.C", "possible.A.to.G", 
									"possible.A.to.T", "possible.A.to.C",
									"possible.T.to.G", "possible.T.to.A", 
									"possible.T.to.C", "possible.C.to.G", 
									"possible.C.to.A", "possible.C.to.T"], 0)
		## Now update sitecounts for this locus.
		## Iterate through each codon in the sequence.
		for codon in codons(dna_sequence):
			sitecounts["possible.G.to.A"] = sitecounts["possible.G.to.A"] + mutate_codon(codon, 'G','A', trans_table)
			sitecounts["possible.G.to.T"] = sitecounts["possible.G.to.T"] + mutate_codon(codon, 'G','T', trans_table)
			sitecounts["possible.G.to.C"] = sitecounts["possible.G.to.C"] + mutate_codon(codon, 'G','C', trans_table)
			sitecounts["possible.A.to.G"] = sitecounts["possible.A.to.G"] + mutate_codon(codon, 'A','G', trans_table)
			sitecounts["possible.A.to.T"] = sitecounts["possible.A.to.T"] + mutate_codon(codon, 'A','T', trans_table)
			sitecounts["possible.A.to.C"] = sitecounts["possible.A.to.C"] + mutate_codon(codon, 'A','C', trans_table)
			sitecounts["possible.T.to.G"] = sitecounts["possible.T.to.G"] + mutate_codon(codon, 'T','G', trans_table)
			sitecounts["possible.T.to.A"] = sitecounts["possible.T.to.A"] + mutate_codon(codon, 'T','A', trans_table)
			sitecounts["possible.T.to.C"] = sitecounts["possible.T.to.C"] + mutate_codon(codon, 'T','C', trans_table)
			sitecounts["possible.C.to.G"] = sitecounts["possible.C.to.G"] + mutate_codon(codon, 'C','G', trans_table)
			sitecounts["possible.C.to.A"] = sitecounts["possible.C.to.A"] + mutate_codon(codon, 'C','A', trans_table)
			sitecounts["possible.C.to.T"] = sitecounts["possible.C.to.T"] + mutate_codon(codon, 'C','T', trans_table)
		## This entry is to check that the number of possible mutations
		## correlates well with gene length (the naive target size.)
		sitecounts["gene.length"] = len(dna_sequence)
		sitecounts["gene"] = gene_name
		locus_dict[locus_tag] = sitecounts

	## The name of the output file is hardcoded here.
	output = writer(open("possible_synonymous.csv",'wb'))
	header = ["locus_tag", "possible.G.to.A", 
			  "possible.G.to.T", "possible.G.to.C",
			  "possible.A.to.G", "possible.A.to.T", "possible.A.to.C",
			  "possible.T.to.G", "possible.T.to.A", "possible.T.to.C",
			  "possible.C.to.G", "possible.C.to.A", "possible.C.to.T",
			  "gene.length", "gene"]
	print header
	output.writerow(header)
	for locus in sorted(locus_dict.keys()):
		cur_dict = locus_dict[locus]
		line_data = [locus, cur_dict["possible.G.to.A"], cur_dict["possible.G.to.T"], cur_dict["possible.G.to.C"],
			   cur_dict["possible.A.to.G"], cur_dict["possible.A.to.T"], cur_dict["possible.A.to.C"],
			   cur_dict["possible.T.to.G"], cur_dict["possible.T.to.A"], cur_dict["possible.T.to.C"],
			   cur_dict["possible.C.to.G"], cur_dict["possible.C.to.A"], cur_dict["possible.C.to.T"], cur_dict["gene.length"], cur_dict["gene"]]
		string_data = [str(x) for x in line_data]
		print string_data
		output.writerow(string_data)
		
main()
