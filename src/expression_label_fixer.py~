#!/usr/bin/python

## expression_label_fixer.py by Rohan Maddamsetti.
## This script fixes the problem that gene names in Tim Cooper's arrays.txt 
## (see his 2003 PNAS paper) don't always match up to the names in my
## ltee_synonymous_mutations.csv.

## This script prints out a list of all genes that map one-to-one between
## REL606.gbk and renamed_arrays.txt (the gene expression data).
## It also print out the total number of synonymous mutations across the LTEE
## by 40,000 generations at these loci.

## Usage: python expression_label_fixer.py > ../data/expression_mutation_data.csv

from Bio import SeqIO


## From stackoverflow, for how to find duplicate elements in a list.
## Useful for finding different loci with the same gene name.
def find_duplicate_elements(seq):
	seen = set()
	seen_add = seen.add
	# adds all elements it doesn't know yet to seen and all other to seen_twice
	seen_twice = set( x for x in seq if x in seen or seen_add(x) )
	# turn the set into a list (as requested)
	return list( seen_twice )

## First, make a list of triples of locus_tags, gene names, and expression names.
## Avoid duplicate entries by hashing by locus_tag.
triple_hash = {}

for genome in SeqIO.parse("/Users/Rohandinho/Desktop/REL606.gbk", "genbank"):
	for feature in genome.features:
		if feature.type != "CDS":
			continue
		else:
			try:
				raw_note = feature.qualifiers["note"][0]
				if raw_note.startswith("b"):
					note = raw_note[0:5] ## disallow actual notes or weirdness.
				else:
					note = ""
			except KeyError:
				note = ""
			try:
				gene = feature.qualifiers["gene"][0]
			except KeyError:
				gene = ""
			locus_tag = feature.qualifiers["locus_tag"][0]
			triple_hash[locus_tag] = (locus_tag, gene, note)
			##print locus_tag, gene, note
            
expression_file = open("/Users/Rohandinho/Desktop/Projects/luscombe-reanalysis/data/renamed_arrays.txt", "r")

## Every entry in database is a gene that maps one-to-one between REL606.gbk
## and the gene expression data.
database = {}

## get the name of each gene in Tim Cooper's gene expression data.
## look up this gene in triple_hash.
for i, line in enumerate(expression_file):
	if i == 0: ## skip the header
		continue
	else:
		name = line.split("\t")[0]
		for locus, triple in triple_hash.iteritems():
			for n in triple:
				if name == n:
					database[locus] = name

## Problem: there are 3548 locus_tags, mapping to 3543 genes. not one-to-one!
## Remove entries in database where the same name is used for multiple locus_tags.
duplicate_genes = find_duplicate_elements(database.values())
bad_entries = []
for k,v in database.iteritems():
	if v in duplicate_genes:
		bad_entries.append(k)
## Now remove the problematic entries.
for i in bad_entries:
	database.pop(i,None)

mutations_file = open("/Users/Rohandinho/Desktop/Projects/luscombe-reanalysis/data/ltee_synonymous_mutations.csv", "r")

## initialize the hash counting synonymous mutations per locus.
mutation_counts = {}
for k in database.keys():
	mutation_counts[k] = 0

for i, line in enumerate(mutations_file):
	if i == 0: ## skip the header
		continue
	else:
		data = line.split(",")
		locus = data[0]
		try:
			mutation_counts[locus] = mutation_counts[locus] + 1
		except KeyError: ## mutation occurred in a gene without expression data.
			continue

print "locus_tag,Gene,synonymous.mutations"
for key in sorted(database.keys()):
	print ",".join([key,database[key], str(mutation_counts[key])])
