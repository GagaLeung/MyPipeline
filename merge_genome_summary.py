#!/usr/bin/env python
"""
@author: Paul W. Bible
@date: 2016-05-14
A tools to merge genome summary data
"""
import sys
import os
import getopt
from collections import defaultdict

clade_order = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
clade_set = set(clade_order)
clade_mark = {"superkingdom":"sk", "kingdom":"k", "phylum":"p", "class":"c", "order":"o",
              "family":"f","genus":"g", "species":"s"}
clade_pred = {"superkingdom":None, "kingdom":"superkingdom", "phylum":"kingdom", "class":"phylum", "order":"class",
              "family":"order","genus":"family", "species":"genus"}


def print_usage(msg=""):
    if len(msg) != 0:
        print msg
    print "usage: merge_genome_summary.py [options] -i <input_directory>"
    print "  looks for a species id in the final column of table and replaces it with the clade path"
    print "  -i, --input=<input_directory>:   A directory containing the output of genome summary analysis."
    print "  -f, --full:                      Print the Full name of the taxon, Default: Only the most specific name"
    print "                                   e.g. Bacteria| ... |Enterococcus|Enterococcus faecalis"
    print "  -t, --taxon_number:              Add the taxon before the species name. e.g. 1351|Enterococcus faecalis"
    print "  -m, --mapped:                    Report the number of sequences mapped. Default: normalized reads"
    print "  -c, --comma                      Use comma separated format. Default tab \\t"


def parse_options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:fmtc', ['help', 'input=', 'full', 'mapped',
                                                             'taxon_number', 'comma'])
    except getopt.GetoptError as err:
        print str(err)
        print_usage()
        sys.exit(2)

    # print usage if not arguments
    if len(opts) == 0 and len(args) == 0:
        print_usage()
        sys.exit(2)

    # initialize default input variables
    input_dir = None
    use_full_name = False
    report_mapped = False
    add_taxon_id = False
    comma_sep = False

    # process arguments
    for opt, arg in opts:
        if opt in ("-h", "-help"):
            print_usage()
            sys.exit(1)
        elif opt in ("-i", "--input"):
            input_dir = arg
        elif opt in ("-f", "--full"):
            use_full_name = True
        elif opt in ("-m", "--mapped"):
            report_mapped = True
        elif opt in ("-t", "--taxon_number"):
            add_taxon_id = True
        elif opt in ("-c", "--comma"):
            comma_sep = True
        else:
            print_usage("Error: Option %s non recognized" % opt)
            sys.exit(1)

    if not input_dir:
        print_usage("Error: an input file is required")
        sys.exit(1)

    kingdom_set = {'b', 'B', 'e', 'E', 'v', 'V', 'a', 'A'}

    # if not kingdom:
    #     kingdom = 'a'
    # elif kingdom and kingdom not in kingdom_set:
    #     print_usage("Error '%s' not recognized as a kingdom: Choose one of %s" % (kingdom, str(kingdom_set)))
    #     sys.exit(1)
    #
    # if anno_level:
    #     anno_levels = anno_level.split(',')
    #     for val in anno_levels:
    #         if val not in clade_set:
    #             print_usage("Error taxonomy level '%s' not recognized: The taxonomy level must be one of %s"
    #                         % (val, str(clade_set)))
    #             sys.exit(1)
    # else:
    #     anno_level = 'all'

    return input_dir, use_full_name, report_mapped, add_taxon_id, comma_sep


def is_numberic(s):
    """
    Take from http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-float-in-python
    :param s:
    :return:
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_format(fname):
    info_map = {"is_good": True, "msg": ""}
    with open(fname) as f:
        header = f.readline()
        for i in xrange(5):
            line = f.readline()
            if not line:
                continue
            row = line.strip().split('\t')
            if len(row) != 5:
                info_map['is_good'] = False
                info_map['msg'] = 'File "%s" does not appear to be a genome summary file' % os.path.basename(fname)
                info_map['msg'] += ' (Number of columns is not 5)'
                break
            elif not is_numberic(row[0]):
                info_map['is_good'] = False
                info_map['msg'] = 'File "%s" does not appear to be a genome summary file' % os.path.basename(fname)
                info_map['msg'] += ' (Number of columns 1 should be digits only)'
                break
            elif not is_numberic(row[2]):
                info_map['is_good'] = False
                info_map['msg'] = 'File "%s" does not appear to be a genome summary file' % os.path.basename(fname)
                info_map['msg'] += ' (Number of columns 3 should be a number)'
                break
            elif not is_numberic(row[3]):
                info_map['is_good'] = False
                info_map['msg'] = 'File "%s" does not appear to be a genome summary file' % os.path.basename(fname)
                info_map['msg'] += ' (Number of columns 4 should be a number)'
                break
            elif not is_numberic(row[4]):
                info_map['is_good'] = False
                info_map['msg'] = 'File "%s" does not appear to be a genome summary file' % os.path.basename(fname)
                info_map['msg'] += ' (Number of columns 4 should be a number)'
                break
    return info_map


def make_in_dir(adir):
    return lambda f: os.path.abspath(adir) + os.path.sep + f


def get_tax_name(use_full_name, tax_name):
    if use_full_name:
        return tax_name
    else:
        return tax_name.split('|')[-1]


def main():
    input_dir, use_full_name, report_mapped, add_taxon_id, comma_sep = parse_options()
    in_data_dir = make_in_dir(input_dir)
    files = os.listdir(input_dir)

    if comma_sep:
        sep_char = ','
    else:
        sep_char = '\t'

    def is_valid_file(a_file):
        format_res = check_format(in_data_dir(a_file))
        if not format_res['is_good']:
            print >> sys.stderr, format_res['msg'], '... skipping'
            return False
        else:
            return True

    good_files = filter(is_valid_file, files)

    norm_read_map = {}
    id_set = set()
    id_to_name = {}
    for a_file in good_files:
        tax_to_reads = {}
        with open(in_data_dir(a_file)) as f:
            header = f.readline()
            count = 0
            for line in f:
                row = line.strip().split('\t')
                tax_id = row[0]
                id_to_name[tax_id] = get_tax_name(use_full_name, row[1])
                id_set.add(row[0])

                if report_mapped:
                    tax_to_reads[tax_id] = row[-2]
                else:
                    tax_to_reads[tax_id] = row[-1]
                count += 1
        norm_read_map[a_file] = tax_to_reads

    id_list = list(id_set)
    id_list.sort(key=lambda x: int(x))

    header = ['species']
    file_name_list = norm_read_map.keys()
    file_name_list.sort()
    for fname in file_name_list:
        header.append(fname)
    print sep_char.join(header)

    # for each taxon
    for i in range(len(id_list)):
        tax_id = id_list[i]
        taxon_label = id_to_name[tax_id]
        if add_taxon_id:
            taxon_label = tax_id + '|' + taxon_label
        row = [taxon_label]
        for fname in file_name_list:
            tax_val_map = norm_read_map[fname]
            if tax_id in tax_val_map:
                row .append(tax_val_map[tax_id])
            else:
                row.append('0')
        print sep_char.join(row)

if __name__ == "__main__":
    main()

