# Usage:

# cd /scratch/steep/7callers/vardict_results

#python ./scripts/vardict_fpfilter.py \
#S1.vcf \
#S1_fpfiltered.vcf

######################################
import sys 

infile=sys.argv[1]

outfile=open(sys.argv[2], 'w')

for line in open(infile):
  if line[0] == '#':
    outfile.write(line)
  if not line[0] == '#':
    columns = line.split('\t')
    status = (columns[7].split(';'))[0].split('=')[1] # Takes the "STATUS=StrongSomatic" string 
    status_title = ((columns[7].split(';'))[0].split('=')[0])
    var_type = (columns[7].split(';')[2].split('=')[1]) # Takes the "TYPE=SNV" string
    var_type_title = (columns[7].split(';')[2].split('=')[0])
    filter_vcf = columns[6]
    DP_title = (columns[8].split(':')[1]) 
    DP = int(columns[9].split(':')[1]) 
    VD_title = (columns[8].split(':')[2])
    VD = int(columns[9].split(':')[2])
    MQ_title = (columns[8].split(':')[14])
    MQ = float(columns[9].split(':')[14])
    SBF_title = (columns[8].split(':')[12])
    SBF = float(columns[9].split(':')[12])
    if status_title != 'STATUS':
        print(line)
        print('Error: STATUS out of position')
        sys.exit()
    if var_type_title != 'TYPE':
        print(line)
        print('Error: VAR TYPE out of position')
        sys.exit()
    if DP_title != 'DP':
        print(line)
        print('Error: DP out of position')
        sys.exit()
    if VD_title != 'VD':
        print(line)
        print('Error: VD out of position')
        sys.exit()
    if MQ_title != 'MQ':
        print(line)
        print('Error: MQ out of position')
        sys.exit()
    if SBF_title != 'SBF':
        print(line)
        print('Error: SBF out of position')
        sys.exit()
    if status == 'StrongSomatic' and (var_type == 'Insertion' or var_type == 'Deletion' or var_type == 'Complex' or var_type == 'SNV') and filter_vcf == 'PASS':
        if DP >= 6:
            if VD >= 2:
                if MQ >= 20.0:
                    if SBF >= 0.01:
                        outfile.write(line)
print('Finished')
######################################
