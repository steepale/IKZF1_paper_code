# Usage:

# cd /scratch/steep/7callers/indelocator_results

# python ./scripts/indelocator_fpfilter.py \
# ./data/S1.vcf \
# ./data/S1_fpfiltered.vcf

# For InDels:
######################################
import sys 
import os
import re

infile=open(sys.argv[1])

outfile=open(sys.argv[2], 'w')

for line in infile:
  if line[0] == '#': 
    outfile.write(line)
  if not line[0] == '#':
  	if re.search(r'SOMATIC', line):
  		columns = line.split('\t')
  		N_DP_title = (columns[7].split(';')[1].split('=')[0])
  		N_DP = int(columns[7].split(';')[1].split('=')[1])
  		N_MQ_title = (columns[7].split(';')[3].split('=')[0])
  		N_MQ = float(columns[7].split(';')[3].split('=')[1].split(',')[1])
  		N_NQSBQ_title = (columns[7].split(';')[4].split('=')[0])
  		N_NQSBQ = float(columns[7].split(';')[4].split('=')[1].split(',')[1])
  		T_AC_title = (columns[7].split(';')[8].split('=')[0])
  		T_AC = int(columns[7].split(';')[8].split('=')[1].split(',')[0])
  		T_DP_title = (columns[7].split(';')[9].split('=')[0])
  		T_DP = int(columns[7].split(';')[9].split('=')[1])
  		T_MQ_title = (columns[7].split(';')[11].split('=')[0])
  		T_MQ = float(columns[7].split(';')[11].split('=')[1].split(',')[0])
  		T_NQSBQ_title = (columns[7].split(';')[12].split('=')[0])
  		T_NQSBQ = float(columns[7].split(';')[12].split('=')[1].split(',')[0])
  		if N_DP_title != 'N_DP':
  			print(line)
  			print('Error: N_DP out of position')
  			sys.exit()
  		if N_MQ_title != 'N_MQ':
  			print(line)
  			print('Error: N_MQ out of position')
  			sys.exit()
  		if N_NQSBQ_title != 'N_NQSBQ':
  			print(line)
  			print('Error: N_NQSBQ out of position')
  			sys.exit()
  		if T_AC_title != 'T_AC':
  			print(line)
  			print('Error: T_AC out of position')
  			sys.exit()
  		if T_DP_title != 'T_DP':
  			print(line)
  			print('Error: T_DP out of position')
  			sys.exit()
  		if T_MQ_title != 'T_MQ':
  			print(line)
  			print('Error: T_MQ out of position')
  		if T_NQSBQ_title != 'T_NQSBQ':
  			print(line)
  			print('Error: T_NQSBQ out of position')
  		if N_DP >= 6:
  			if N_MQ >= 20.0:
  				if N_NQSBQ >= 25.0:
  					if T_AC >= 2:
  						if T_DP >= 6:
  							if T_MQ >= 20.0:
  								if T_NQSBQ >= 25.0:
  									#print(pvalue)
  									outfile.write(line)
print('Finished')
######################################




