#!/usr/bin/python

## tabulateSynonymousInLTEE.py by Rohan Maddamsetti.
## This script takes a set of annotated genomediff files,
## and it produces a csv file tabulating locus_tag, the genome source,
## and the number of G.to.A (etc.) mutations at that locus in the given
## genome.
## Usage: python tabulateSynonymousInLTEE.py ../data/40K_annotated_diffs/*.gd

from sys import argv
from os.path import basename
from csv import writer
import sys

def main():
	list_of_input_gds = argv[1:]
	genome_dict = {}
	for gd_file in list_of_input_gds:
		##get the name of the genome.
		genome_name = basename(gd_file).split('.')[0]
		gd_handle = open(gd_file)
		genome_dict[genome_name] = {}
		for line in gd_handle:
			if "snp_type=synonymous" in line:
				data = line.split("\t")
				gene_name = data[13].split('=')[1]
				locus_tag = data[17].split('=')[1]
				new_codon = data[9].split('=')[1]
				old_codon = data[11].split('=')[1]
				##figure out what the new mutation is, from the codons.
				old,new = None, None
				for i,j in zip(old_codon,new_codon):
					if i != j:
						old, new = i,j
						break
				##save these data in some data structure.
				gene_dict = genome_dict[genome_name]
				if locus_tag not in gene_dict:
					##initialize the entry (without the gene name entry yet)
					gene_dict[locus_tag] = dict.fromkeys(["G.to.A", "G.to.T", 
														  "G.to.C", "A.to.G", 
														  "A.to.T", "A.to.C",
														  "T.to.G", "T.to.A", 
														  "T.to.C", "C.to.G", 
														  "C.to.A", "C.to.T"], 0)
				## Now update the entry.
				sitecounts = gene_dict[locus_tag]
				sitecounts["gene_name"] = gene_name
				if old == 'G':
					if new == 'A':
						sitecounts["G.to.A"] = sitecounts["G.to.A"] + 1
					if new == 'T':
						sitecounts["G.to.T"] = sitecounts["G.to.T"] + 1
					if new == 'C':
						sitecounts["G.to.C"] = sitecounts["G.to.C"] + 1
				elif old == 'A':
					if new == 'G':
						sitecounts["A.to.G"] = sitecounts["A.to.G"] + 1
					if new == 'T':
						sitecounts["A.to.T"] = sitecounts["A.to.T"] + 1
					if new == 'C':
						sitecounts["A.to.C"] = sitecounts["A.to.C"] + 1
				elif old == 'T':
					if new == 'G':
						sitecounts["T.to.G"] = sitecounts["T.to.G"] + 1
					if new == 'A':
						sitecounts["T.to.A"] = sitecounts["T.to.A"] + 1
					if new == 'C':
						sitecounts["T.to.C"] = sitecounts["T.to.C"] + 1
				elif old == 'C':
					if new == 'G':
						sitecounts["C.to.G"] = sitecounts["C.to.G"] + 1
					if new == 'A':
						sitecounts["C.to.A"] = sitecounts["C.to.A"] + 1
					if new == 'T':
						sitecounts["C.to.T"] = sitecounts["C.to.T"] + 1
			else:
				continue

	## Dump the data into a csv file.
	output = writer(open("ltee_synonymous_mutations.csv","wb"))
	header = ["locus_tag", "gene", "G.to.A", "G.to.T", "G.to.C",
			  "A.to.G", "A.to.T", "A.to.C", 
			  "T.to.G", "T.to.A", "T.to.C",
			  "C.to.G", "C.to.A", "C.to.T", "genome"]
	output.writerow(header)
	for this_genome, this_gene_dict in genome_dict.iteritems():
		for this_locus, fields in this_gene_dict.iteritems():
			rowdata = [this_locus, fields["gene_name"], fields["G.to.A"], fields["G.to.T"], 
					  fields["G.to.C"], fields["A.to.G"], fields["A.to.T"], 
					  fields["A.to.C"], fields["T.to.G"], fields["T.to.A"],
					  fields["T.to.C"],	fields["C.to.G"], fields["C.to.A"],
					  fields["C.to.T"], this_genome]
			output.writerow(rowdata)


main()
