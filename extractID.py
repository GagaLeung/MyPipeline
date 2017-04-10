# import os
# os.chdir('/Users/gaga/Desktop/filter_sample')
# bed = 'ERR011103.coverage.txt'
import sys
bed = sys.argv[1]  # output txt from genomecovBed
cov = open(bed)
# column 1 : scaffold number
# column 2 : covered times
# column 3 : covered length
# column 4 : reference length
# column 5 : scaffold coverage (column 3 / column 4)
output_name = 'extracted_' + bed[0:9] + '.txt'
ext = open(output_name, 'w')
for line in cov:
    c1 = line.split('\t')[0]
    if c1 == 'genome':
        continue
    c2 = line.split('\t')[1]
    if int(c2) != 0:
        continue
    c3 = line.split('\t')[2]
    c4 = line.split('\t')[3]
    scaf = c1.split('|')
    taxid = scaf[1]
    new = taxid + '\t' + c3 + '\t' + c4 + '\n'
    ext.write(new)

ext.close()

