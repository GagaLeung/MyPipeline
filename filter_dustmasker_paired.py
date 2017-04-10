#!/bin/env python
"""
@author: Paul W. Bible
@date: 2016-10-10
A tool for filtering out sequences with low complexity from dustmasker
"""
import sys
import getopt


def next_record(handle):
    line = handle.readline()
    # print "next_fasta call:", line
    if not line or line[0] != ">":
        return None

    header = line.strip('>').split()[0]
    rec_buffer = []

    current = handle.tell()
    line = handle.readline()
    while len(line) != 0 and line[0] != ">":
        rec_buffer.append(line.strip())
        current = handle.tell()
        line = handle.readline()
    handle.seek(current)

    return header, rec_buffer


def rec_sequence(handle):
    rec = next_record(handle)
    while rec:
        yield rec
        rec = next_record(handle)


def shared_prefix(filename_1, filename_2):
    index = 0
    for i in range(len(filename_1)):
        if filename_1[i] != filename_2[i]:
            index = i - 1
            break
    base = filename_1[:index]
    base = base.strip('_')
    base = base.strip('.')
    return base


def print_usage(msg=""):
    if len(msg) != 0:
        print msg
    print 'usage: filter_dustmasker_paired.py -1 <input_1> -2 <input_2> [-l <limit>]'
    print '  Create process scripts for a directory of paired end data'
    print
    print ' -1                  The first dustmasker input file.'
    print ' -2                  The second dustmasker input file.'
    print ' -l <limit>          Mask limit, filter sequences with greater than <limit> bases, default is 30'


def process_and_validate_arguments():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h1:2:', ['help','input1=', 'input2=', 'limit='])
        #print opts
        #print args
    except getopt.GetoptError as err:
        print str(err)
        print_usage()
        sys.exit(2)

    if len(opts) == 0:
        print_usage()
        sys.exit(1)

    input1 = None
    input2 = None
    limit = None

    for opt, arg in opts:
        if opt in ("-h","-help", "-help"):
            print_usage()
            sys.exit(1)
        elif opt in ("-1", "--input1"):
            input1 = arg
        elif opt in ("-2", "--input2"):
            input2 = arg
        elif opt in ("-l", "--limit"):
            limit = arg
        else:
            print_usage("Error: Option %s non recognized" % opt)
            sys.exit(1)

    if not input1:
        print_usage("Error: Missing input 1, use '-1 <input_1>'")
        sys.exit(2)

    if not input2:
        print_usage("Error: Missing input 2, use '-2 <input_2>'")
        sys.exit(2)

    if limit and not limit.isdigit():
        print_usage("Error: Limit must be a digit.")
        sys.exit(2)

    return input1, input2, limit


def xor(val1, val2):
    return bool(val1) != bool(val2)


def calculate_total_masked_bases(masks):
    total = 0
    for mask in masks:
        start, end = mask.split('-')
        total += int(end) - int(start)
    return total


def main():
    input1, input2, limit = process_and_validate_arguments()
    if limit:
        threshold = int(limit)
    else:
        threshold = 30

    with open(input1, 'rb') as f, open(input2, 'rb') as f2:
        for header, masks in rec_sequence(f):
            header2, masks2 = next_record(f2)
            if len(masks) > 0 and len(masks2) > 0:
                total1 = calculate_total_masked_bases(masks)
                total2 = calculate_total_masked_bases(masks)
                if total1 < threshold and total2 < threshold:
                    print header
                    print >> sys.stderr, header2
            elif xor(len(masks) > 0, len(masks2) > 0):
                #print header, header2, masks, masks2
                combined = masks
                combined.extend(masks2)
                total = calculate_total_masked_bases(combined)
                if total < threshold:
                    print header
                    print >> sys.stderr, header2
            else:
                print header
                print >> sys.stderr, header2


if __name__ == '__main__':
    main()

