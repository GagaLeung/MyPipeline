import sys
sample = sys.argv[1]
# import os
# os.chdir('/Users/gaga/Desktop/filter_sample')
# cut1 = 'cut0_' + sample + '.txt'
# cut2 = 'cut0.01_' + sample + '.txt'
# cut3 = 'cut0.1_' + sample + '.txt'
# cut4 = 'cut1_' + sample + '.txt'
# cut5 = 'cut5_' + sample + '.txt'
# cut6 = 'cut10_' + sample + '.txt'
# cut7 = 'cut20_' + sample + '.txt'


def beautify(access, cutoff):
    # access is sample name from command line;
    # cutoff is a number that can be one of 0, 0.01, 0.1, 1, 5, 10, 20
    r_out = 'cut' + str(cutoff) + '_' + access + '.txt'
    fh_from = open(r_out)
    py_out = 'all_format_' + r_out
    fh_to = open(py_out, 'w')
    fh_to.write('taxon_id' + '\t' + 'taxonomy_names' + '\t' + 'avg_genome_length' + '\t' + 'reads_mapped' + '\t' + 'length_normalized_reads' + '\n')
    for line in fh_from:
        if line.startswith('"taxid"'):
            continue
        first = line.find(' ')
        second = line.find(' ', first + 1)
        third = line.find(' ', second + 1)
        taxid = line[first + 1: second]
        if int(taxid) != 12355 and int(taxid) != 322067:
            quo_after_tax_name = line.find('"', third)
        else:
            quo_after_tax_name = line.find('"', second + 2)
        blk_after_gl = line.find(' ', quo_after_tax_name + 2)
        blk_after_rm = line.find(' ', blk_after_gl + 1)
        blk_after_ln = line.find(' ', blk_after_rm + 1)

        tax_name = line[second + 2: quo_after_tax_name]
        avg_genome_length = line[quo_after_tax_name + 2: blk_after_gl]
        reads_mapped = line[blk_after_gl + 1: blk_after_rm]
        length_normalized_reads = line[blk_after_rm + 1: blk_after_ln]
        fh_to.write(
                taxid + '\t' + tax_name + '\t' + avg_genome_length + '\t' + reads_mapped + '\t' + length_normalized_reads + '\n')

for cut in [0, 0.01, 0.1, 1, 5, 10, 20]:
    beautify(sample, cut)




