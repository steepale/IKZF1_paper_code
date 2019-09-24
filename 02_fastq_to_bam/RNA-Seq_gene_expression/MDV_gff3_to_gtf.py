import os
import sys

infile = sys.argv[1]
outfile = open(sys.argv[2], 'w')
gene_table = open(sys.argv[3], 'w')

# The gff3 and gtf files are conveniently ordered by feature as 
# gene, transcript, exon, cds or some subset but always in that hierarchy

# Variables
keeps = ('gene','transcript','exon','CDS')

# Write a header for the gene table file
gene_table.write('\t'.join(['GENE_ID','PRODUCT'])+'\n')
gene_table.write('\t'.join(['5291091','vTR(5291091):telomerase_RNA_subunit'])+'\n')
gene_table.write('\t'.join(['5291094','vTR(5291094):telomerase_RNA_subunit'])+'\n')

# Make one pass and construct a dictionary
# The gene ID's or parent ID's will be the keys, and the values will be lists
for line in open(infile):
	if line.startswith("##"):
		outfile.write(line)
	if not line.startswith("#"):
		col = line.rstrip().split('\t')
		CHROM = col[0]
		ANNO = col[1]
		FEATURE = col[2]
		START = col[3]
		END = col[4]
		SCORE = col[5]
		STRAND = col[6]
		FRAME = col[7]
		ATTR = col[8]
		# Only collect attributes for genes and features associated with genes
		if FEATURE in keeps:
			if FEATURE == 'gene':
				# When a new gene is encountered, reset all variables
				gene_id,gene_source,gene_biotype,gene_name,transcript_id,transcript_source,transcript_biotype,exon_id,protein_id = ["NA"] * 9
				# Assign: gene_id, gene_biotype, gene_name, gene_source
				gene_id = ATTR.split('Dbxref=GeneID:')[1].split(';')[0] # Absolute
				gene_source = ANNO # Absolute
				gene_biotype = ATTR.split('gene_biotype=')[1].split(';')[0] # Absolute
				gene_name = ATTR.split('gene=')[1].split(';')[0] # Conditional
				Name = ATTR.split('Name=')[1].split(';')[0]
				if gene_name != "NA": # If gene_name is assigned
					GENE_ATTR = 'gene_id "'+gene_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'gene_name "'+gene_name+'";'
				elif gene_name == "NA": # If gene_name is not assigned
					GENE_ATTR = 'gene_id "'+gene_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'";'
				outfile.write('\t'.join([CHROM,ANNO,FEATURE,START,END,SCORE,STRAND,FRAME,GENE_ATTR])+'\n')
			elif FEATURE == 'transcript':
				# Assign: transcript_id, Parent, transcript_source
				transcript_id = ATTR.split('ID=')[1].split(';')[0]
				transcript_source = ANNO
				transcript_biotype = gene_biotype # There is no transcript biotype annotation in MDV gff3 file, so we will use gene
				# Absolutes: gene_id (1), transcript_id(2), gene_source(3), gene_biotype(4),transcript_source(5), transcript_biotype(6)
				# Conditionals: gene_name(2.1)
				if gene_name != "NA": # If gene_name is assigned
					TRANS_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_name "'+gene_name+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'";'
				elif gene_name == "NA": # If gene_name is not assigned
					TRANS_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'";'
				outfile.write('\t'.join([CHROM,ANNO,FEATURE,START,END,SCORE,STRAND,FRAME,TRANS_ATTR])+'\n')
			elif FEATURE == 'exon':
				exon_id = ATTR.split('ID=')[1].split(';')[0]
				# Absolutes: gene_id(1), transcript_id(2), gene_source(4), gene_biotype(5), transcript_source(6), transcript_biotype(7),exon_id(8)
				# Conditionals: gene_name(3)
				if gene_name != "NA": # If gene_name is assigned
					EXON_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_name "'+gene_name+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'exon_id "'+exon_id+'";'
				elif gene_name == "NA": # If gene_name is not assigned
					EXON_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'exon_id "'+exon_id+'";'
				outfile.write('\t'.join([CHROM,ANNO,FEATURE,START,END,SCORE,STRAND,FRAME,EXON_ATTR])+'\n')
			elif FEATURE == 'CDS':
				# Absolutes: gene_id(1), transcript_id(2), gene_source(4), gene_biotype(5), transcript_source(6), transcript_biotype(7),protein_id(8),
				# Conditionals: gene_name(3)
				protein_id = ATTR.split('protein_id=')[1].split(';')[0]
				product = ATTR.split('product=')[1].split(';')[0].replace(' ','_')
				product = Name+'('+gene_id+'):'+product
				# HTSEQ will only count 1 feature type. The default type is exon and is prefered for ensemble gtf files.
				# We will merge this file with the chicken ensembl file and therefore we will adjust each CDS to also be considered as an exon for counting purposes.
				if gene_name != "NA": # If gene_name is assigned
					EXON_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_name "'+gene_name+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'exon_id "'+exon_id+'";'
				elif gene_name == "NA": # If gene_name is not assigned
					EXON_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'exon_id "'+exon_id+'";'
				outfile.write('\t'.join([CHROM,ANNO,'exon',START,END,SCORE,STRAND,FRAME,EXON_ATTR])+'\n')
				if gene_name != "NA": # If gene_name is assigned
					CDS_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_name "'+gene_name+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'protein_id "'+protein_id+'";'
				elif gene_name == "NA": # If gene_name is not assigned
					CDS_ATTR = 'gene_id "'+gene_id+'"; '+'transcript_id "'+transcript_id+'"; '+'gene_source "'+gene_source+'"; '+'gene_biotype "'+gene_biotype+'"; '+'transcript_source "'+transcript_source+'"; '+'transcript_biotype "'+transcript_biotype+'"; '+'protein_id "'+protein_id+'";'
				outfile.write('\t'.join([CHROM,ANNO,FEATURE,START,END,SCORE,STRAND,FRAME,CDS_ATTR])+'\n')
				# Create another file with gene_id and gene description
				gene_table.write('\t'.join([gene_id,product])+'\n')
gene_table.close()
outfile.close()
