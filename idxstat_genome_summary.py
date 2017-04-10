#!/usr/bin/env python
"""
@author: Paul W. Bible
@date: 2016-05-14
A tool for extracting the full clade path of a speceis
"""
import sys
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
    print "usage: species_path.py [options] -n <node_file> -m <name_file> -i <input_file>"
    print "  looks for a species id in the final column of table and replaces it with the clade path"
    print "  -n, --nodes=<node_file>:    The taxonomy file for NODES (nodes.dmp) from NCBI"
    print "  -m, --names=<names_file>:   The taxonomy file for NAMES (names.dmp) from NCBI"
    print "  -i, --input=<input_file>:   an input table (tab delimited) with the last column a species or taxon id."
    print "  -l, --level:                The annotation level. Print only annotations at this level"
    print "                              Levels can be one of kingdon, phylum, class, family, genus, or species"
    print "  -k, --kingdom:              Limit analysis to this kingdom only, values are:"
    print "                                  'b' for Bacteria, the default"
    print "                                  'e' for Eukaryota"
    print "                                  'v' for Virus"
    print "                                  'a' for All kingdoms, calculates data for all kingdoms, maybe skewed by" \
          "small genomes."


def parse_options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hk:l:n:m:i:a', ['help', 'nodes=', 'names=', 'input=',
                                                                'annotate', 'level', 'kingdom'])
    except getopt.GetoptError as err:
        print str(err)
        print_usage()
        sys.exit(2)

    if len(opts) == 0 and len(args) == 0:
        print_usage()
        sys.exit(2)

    nodes_file = ""
    names_file = ""
    input_file = ""
    anno_level = None
    annotate = False
    kingdom = None

    for opt, arg in opts:
        if opt in ("-h","-help"):
            print_usage()
            sys.exit(1)
        elif opt in ("-n", "--nodes"):
            nodes_file = arg
        elif opt in ("-m", "--names"):
            names_file = arg
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-a", "--annotate"):
            annotate = True
        elif opt in ("-l", "--level"):
            anno_level = arg
        elif opt in ("-k", "--kingdom"):
            kingdom = arg
        else:
            print_usage("Error: Option %s non recognized" % opt)
            sys.exit(1)

    if nodes_file == "":
        print_usage("Error: nodes file is required")
        sys.exit(1)

    if names_file == "":
        print_usage("Error: names file is required")
        sys.exit(1)

    if input_file == "":
        print_usage("Error: an input file is required")
        sys.exit(1)

    kingdom_set = {'b', 'B', 'e', 'E', 'v', 'V', 'a', 'A'}

    if not kingdom:
        kingdom = 'a'
    elif kingdom and kingdom not in kingdom_set:
        print_usage("Error '%s' not recognized as a kingdom: Choose one of %s" % (kingdom, str(kingdom_set)))
        sys.exit(1)

    if anno_level:
        anno_levels = anno_level.split(',')
        for val in anno_levels:
            if val not in clade_set:
                print_usage("Error taxonomy level '%s' not recognized: The taxonomy level must be one of %s"
                            % (val, str(clade_set)))
                sys.exit(1)
    else:
        anno_level = 'all'

    return nodes_file, names_file, input_file, annotate, anno_level, kingdom


def parse_nodes_file(fname):
    pred_map = {}
    clade_type_map = {}
    f = open(fname)
    for line in f:
        line = line.strip()
        # print line
        items = line.split("\t|\t")
        # print items
        # print len(items)
        node_id = int(items[0])
        parent_id = int(items[1])
        clade_type = items[2]
        pred_map[node_id] = parent_id
        clade_type_map[node_id] = clade_type
    f.close()
    return pred_map, clade_type_map


def parse_names_file(fname):
    name_map = {}
    f = open(fname)
    for line in f:
        line = line.strip()
        if line.find("scientific name") < 0:
            continue

        line = line.rstrip("\t|")
        items = line.split("\t|\t")
        node_id = int(items[0])
        clade_name = items[1]
        name_map[node_id] = clade_name
    f.close()
    return name_map


def species_to_clade_map(species_id, parent_map, clade_levels_map, name_map):
    levels = []
    names = []

    curr_id = species_id
    species_tax_level_map = {}
    while curr_id != 1:
        if curr_id not in clade_levels_map:
            return "no_path_" + str(curr_id)
        # print curr_id
        clade_level = clade_levels_map[curr_id]
        if clade_level == "no rank" or clade_level not in clade_set:
            curr_id = parent_map[curr_id]
            continue
        # print clade_level
        levels.append(clade_level)

        level_name = name_map[curr_id]
        names.append(level_name)
        species_tax_level_map[clade_level] = level_name

        # print level_name
        curr_id = parent_map[curr_id]
    species_tax_level_map['tax'] = species_id
    return species_tax_level_map


def best_tax_level(level, s_map):
    '''
    Some taxonomy levels are missing important classification nodes. This function will try to fill them with the
    next highest level.
    :param level:
    :param s_map:
    :return:
    '''
    current_level = level
    while current_level not in s_map and current_level != 'superkingdom':
        current_level = clade_pred[current_level]
    if current_level in s_map:
        return s_map[current_level]
    else:
        return 'unknown'


def order_species_map(s_map):
    names = []
    for val in clade_order:
        if val not in s_map:
            names.append(best_tax_level(val, s_map))
        else:
            names.append(s_map[val])
    return names


def genome_id(rec):
    parts = rec[0].split('|')
    return parts[-1]


def tax_id(rec):
    parts = rec[0].split('|')
    return parts[-2]


def seq_length(rec):
    return int(rec[1])


def map_count(rec):
    return int(rec[2])


def main():
    nodes_fname, names_fname, input_file, do_marks, anno_level, kingdom = parse_options()
    parent_map, clade_levels = parse_nodes_file(nodes_fname)
    name_map = parse_names_file(names_fname)

    # Collect read data for each sequences and combine counts by genome
    genome_map = defaultdict(list)
    genome_count_map = defaultdict(int)
    reads_mapped = 0
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            rec = line.split("\t")
            # skip unmapped reads
            if rec[0] == '*':
                continue
            # Collect genome data
            g_id = genome_id(rec)
            genome_map[g_id].append(rec)
            # get the number of mapped reads to this sequence and add to its genome
            genome_count_map[g_id] += map_count(rec)
            # add mapped reads to the total mapped reads of all sequences
            reads_mapped += map_count(rec)

    # print 'Total Reads Mapped:', reads_mapped

    # sort genomes by read count
    vals = genome_count_map.items()
    vals.sort(key=lambda x: -x[1])

    unique_species_map = defaultdict(list)
    unique_species = set()
    total_cpm = 0
    for k, v in vals:
        if v == 0:
            continue

        seqs = genome_map[k]
        genome_length = sum(map(seq_length, seqs))
        rec = seqs[0]
        tax = int(tax_id(rec))
        # print 'genome', k
        # print 'mapped', v
        # print 'genome length', genome_length
        map_ratio = float(v)/genome_length

        # print 'reads per bp', map_ratio
        # species_path = species_to_clade_path(tax, parent_map, clade_levels, name_map, anno_level, mark=do_marks)
        species_path_map = species_to_clade_map(tax, parent_map, clade_levels, name_map)

        # print order_species_map(species_path_map)
        # raw_input()

        # if species_path.find('Viruses') >= 0:
        #     continue
        # print species_path
        # unique_species.add(species_path)
        cpm = (map_ratio * reads_mapped)/1000000
        total_cpm += cpm
        out_rec = [k, 'taxon_' + str(tax), str(genome_length), str(v), str(cpm)]

        # print '\t'.join(out_rec)
        # analyze_recs.append(out_rec)

        if type(species_path_map) == str:
            unique_species_map[species_path_map].append((out_rec, species_path_map))
        else:
            unique_species_map[species_path_map['species']].append((out_rec, species_path_map))



    # Aggregate data at the species level
    species_count_recs = []
    for s in unique_species_map:
        # print s
        genome_recs = unique_species_map[s]
        total_genome_lengths = 0
        total_genome_reads = 0
        # Collect the genome data by species
        for out_rec, s_map in genome_recs:
            # print out_rec, order_species_map(s_map)
            total_genome_lengths += int(out_rec[2])
            total_genome_reads += int(out_rec[3])
        # Calculate the average genome length
        avg_genome_size = float(total_genome_lengths) / len(genome_recs)
        # Format the desired output record

        if anno_level == 'all':
            names = order_species_map(s_map)
        else:
            parts = anno_level.split(',')
            names = [best_tax_level(val, s_map) for val in parts]
        # names.insert(0, str(s_map['tax']))
        species_ids = '|'.join(names)

        # new_rec = [s_map['species'], str(int(avg_genome_size)), str(total_reads)]
        map_ratio = float(total_genome_reads)/avg_genome_size
        cpm = (map_ratio * reads_mapped)
        if type(s_map) == str:
            taxid = s_map.split('_')[-1]
            new_rec = [taxid, species_ids, str(int(avg_genome_size)), str(total_genome_reads), str(cpm)]
        else:
            new_rec = [str(s_map['tax']), species_ids, str(int(avg_genome_size)), str(total_genome_reads), str(cpm)]

        # Only output the specified kingdom
        if kingdom in {'a', 'A'}:
            species_count_recs.append(new_rec)
        elif kingdom in {'b', 'B'} and best_tax_level('kingdom', s_map) == 'Bacteria':
            species_count_recs.append(new_rec)
        elif kingdom in {'e', 'E'} and best_tax_level('kingdom', s_map) == 'Eukaryota':
            species_count_recs.append(new_rec)
        elif kingdom in {'v', 'V'} and best_tax_level('kingdom', s_map) == 'Viruses':
            species_count_recs.append(new_rec)

    header = ['taxon_id','taxonomy_names', 'avg_genome_length', 'reads_mapped', 'length_normalized_reads']
    species_count_recs.sort(key=lambda x: -float(x[-1]))
    print '\t'.join(header)
    for rec in species_count_recs:
        print '\t'.join(rec)


if __name__ == "__main__":
    main()

