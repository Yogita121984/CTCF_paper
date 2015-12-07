#!/usr/bin/env python
# extract_common.py --- 
# 
# Filename: extract_common.py
# Description: This script extracts the common reads from two files.
#
# Usage: extract_common.py [-h] [-c COMMON_FILE] [-p UNIQUE_PREFIX] FILE1 FILE2
# 
# This script extracts the common reads from two files.
# 
# positional arguments:
#   FILE1             First file to read
#   FILE2             Second file to read
#  
# optional arguments:
#   -h, --help        show this help message and exit
#   -c COMMON_FILE    Output file where common reads are saved
#   -p UNIQUE_PREFIX  Prefix appended to unique reads in each file
# 
# Common file format: (co_ORIG_FILE_NAME)
# 
# <chromosome> <start_site> <end_site> <REST_OF_LINE IN ORIGINAL FILE> 
# 
# Overlap file format (common.txt): 
# 
# <chromosome> <overlapping start> <overlapping end> <FILE1 region start> <FILE1 region end> <FILE2 region start> <FILE2 region end> 
#   
# The regions are tab-separated while within a region, separated by a space.  
#
# Created: Tue Aug 13 09:26:53 2013 (+0200)
# Version: 2.0 
# Last-Updated: Mon Oct 28 07:41:15 2013 (+0100)
#           By: Vinay Jethava
#     Update #: 253
# 
# Change Log:
# 27-Oct-2013    Vinay Jethava  
#    Last-Updated: Sun Oct 27 21:46:36 2013 (+0100) #224 (Vinay Jethava)
#    Generates files with co_ prefix which contain the original regions for each file which overlap.
#    Added rest_of_line to dictionary
#    Added logging 
# 
# 17-Oct-2013    Vinay Jethava  
#    Added removal of empty line from script 
# 

import argparse
import sys
import re 
import numpy

from logger import logger 
from IntervalTree import Interval,IntervalTree

# Argument parser
parser = argparse.ArgumentParser(description='This script extracts the common reads from two files.')
parser.add_argument('f1', metavar='FILE1', type=str, help='First file to read')
parser.add_argument('f2', metavar='FILE2', type=str, help='Second file to read') 
parser.add_argument('-c', metavar='COMMON_PREFIX', type=str, help='Prefix appended to output file for common reads (default: co_)', default='co_')
parser.add_argument('-p', metavar='UNIQUE_PREFIX', type=str, help='Prefix appended to unique reads in file (default: u_)', default='u_') 
parser.add_argument('-o', metavar='OVERLAP_FILE', type=str, help='File containing region overlaps (default: overlap.txt)', default='overlap.txt')
parser.add_argument('-v', metavar='VERBOSITY', help="Increase verbosity level [0, 1, 2] (default: 1)", type=int, default=1) 

# separators and joiners (file formatting)
my_separator = '_' 
my_joiner='\t'

def read_file_to_set(fid, separator=my_separator):
    lines = []
    fid.seek(0) 
    for line in fid.readlines(): 
        line = line.strip()
        parts = re.split('\s' , line)
        line = separator.join(parts[0:3])
        print line 
        lines.append(line) 
    return set(lines) 

def read_file_to_map(fid, separator=my_separator): 
    """ 
    Converts output of read_file_to_set to required dict format 

    string => [ (int, int) ] 
    """
    lines = []
    result = dict()
    fid.seek(0) 
    for line in fid.readlines(): 
        line = line.strip()
        parts = re.split('\s' , line )
        parts = filter(None, parts) 
        key = parts[0] # chromosome name
        try: 
            value = (int(round(float(parts[1]))) , int(round(float(parts[2])))) # (start, end)
            logger.debug("file: %s key: %s value: %s" % (fid.name, key, str(value)))
        except:             
            print >> sys.stderr, "line is: ", line
            for i in range(len(parts)): 
                print >> sys.stderr, "parts[%d]: %s" % (i, parts[i])
            raise
        addendum = '\t'.join(parts[3:])
        logger.debug("file: %s key: %s value: %s" % (fid.name, key, str(value)) )
        logger.debug("addendum: %s" % (addendum) )
        if key in result: 
            result[key].append((value, addendum)) 
        else:
            result[key] = [(value, addendum)] 
    for key in result: 
        logger.info("%s chromosome: %s reads: %d" % (fid.name, key,  len(result[key])) )
    return result

############################################################
## Uses IntervalTree O(n log(m) )
############################################################
def find_overlapping_sequences(list1, list2): 
    m = len(list1) 
    n = len(list2) 
    if m < n:
        (result_swapped, uniq2, uniq1) = find_overlapping_sequences(list2, list1)
        result = dict()
        for item in result_swapped.keys():
            i = item[0]
            j = item[1]
            result[(j, i)] = result_swapped[item]
        return (result, uniq1, uniq2)
    else:
        intervals1 = [ Interval( list1[i][0][0], list1[i][0][1], i) for i in range(m)]
        intervals2 = [ Interval( list2[i][0][0], list2[i][0][1], i) for i in range(n)]
        result = dict()
        unique1 = set(range(m))
        unique2 = set(range(n)) 
        interval_tree = IntervalTree.init_from_list(intervals1)
        for i in range(n):
            item2 = intervals2[i]
            overlapping_intervals = interval_tree.overlapping_interval_search(item2)
            for item1 in overlapping_intervals:
                curr_overlap = item1.find_overlap(item2)                
                result[ (item1.id, item2.id) ] = (curr_overlap.left, curr_overlap.right)
                unique1.discard(item1.id)
                unique2.discard(item2.id)
        return (result, unique1, unique2) 


############################################################
## Deprecated version O(m x n)
############################################################
def find_overlapping_sequences_deprecated(list1, list2): 
    """
    Given pairs of (start, end) of reads, finds overlapping reads. 
    by returning a graph (directed) in form of adjacency matrix where 

    A(i, j) = 1 if time_period_i in list1 and time_period_j in list2 overlap
    
    """    
    m = len(list1) 
    n = len(list2) 
    result = dict()
    unique1 = set(range(m))
    unique2 = set(range(n)) 
    
    for i in range(m): 
         for j in range(n):
             overlap = find_overlap(list1[i], list2[j])
             if not (overlap == None):
                 result[(i, j)] = overlap
                 unique1.discard(i)
                 unique2.discard(j)
    return (result, unique1, unique2)  

def find_overlap(period1, period2):
    b1 = period1[0][0]
    e1 = period1[0][1]
    b2 = period2[0][0]
    e2 = period2[0][1]
    overlap = None
    # if overlap: 
    #     return 1
    # else:
    #     return 0
    if ((e1 > b2) and (b1 < e2)) or ((b1 < e2) and (b2 < e1)): 
        overlap =  (max(b1, b2), min(e1, e2))
    return overlap


    
def write_set_to_file(s, fid, separator=my_separator, joiner=my_joiner): 
    for line in s: 
        parts = re.split(separator, line)
        line = joiner.join(parts) 
        print >> fid, line

def write_unique_sites_per_chromosome_to_file(fid, chromosome_name, list_of_chromosome_sites, unique_ids): 
    logger.info("Writing %d unique regions on chromosome: %s to file: %s" % (len(unique_ids), chromosome_name, fid.name)) 
    for key in unique_ids: 
        start_site = list_of_chromosome_sites[key][0][0]
        end_site = list_of_chromosome_sites[key][0][1]
        rest_of_line = list_of_chromosome_sites[key][1]
        fid.write("%s\t%d\t%d\t%s\n" % (chromosome_name, start_site, end_site, rest_of_line))
        logger.debug("file: %s unique: %s (%s, %s)" % (fid.name, chromosome_name, start_site, end_site)) 

if __name__ == '__main__':
    try:
        args = parser.parse_args(sys.argv[1:]) 
    except: 
        print "\n" 
        # parser.print_help()
        raise

    # set verbosity.
    if args.v < 1:
        logger.setLevel(logging.ERROR)
        channel.setLevel(logging.ERROR) 
    elif args.v > 1:
        logger.setLevel(logging.DEBUG)
        channel.setLevel(logging.DEBUG)

    file1 = args.f1 
    file2 = args.f2 
    prefix = args.p 
    # overlap file 
    common_file = args.o 
    # unique files
    unique1_file = prefix + file1 
    unique2_file = prefix + file2 
    # original lines which are common. 
    common_prefix = args.c
    common1_file = common_prefix + file1
    common2_file = common_prefix + file2 
    
    # Read the first file
    with open(file1, 'r') as fid1: 
        #    f1s = read_file_to_set(fid1)
        f1m = read_file_to_map(fid1) 
        
    # Read the second file
    with open(file2, 'r') as fid2: 
        f2m = read_file_to_map(fid2) 
        #       f2s = read_file_to_set(fid2) 

    chromosomes1 = f1m.keys()
    chromosomes2 = f2m.keys()
    common_chr = set(chromosomes1).intersection(set(chromosomes2))

    try: 
        cfid = open(common_file, 'w')
        u1fid = open(unique1_file, 'w') 
        u2fid = open(unique2_file, 'w') 

        ## Writing full lines for each of the common regions.  
        c1fid = open(common1_file, 'w') 
        c2fid = open(common2_file, 'w') 
        
    except:
        print "Could  not open file for writing: ", common_file 
        raise 
              
    for chr in common_chr:
        list1 = f1m[chr]
        list2 = f2m[chr] 
        (overlap_list, unique1, unique2) = find_overlapping_sequences(list1, list2) 

        ############################################################
        ## Comparison with find_overlapping_sequences_deprecated()
        ############################################################
        # (old_ol, old_u1, old_u2) = find_overlapping_sequences_deprecated(list1, list2)
        # assert(overlap_list == old_ol)
        # assert(unique1 == old_u1)
        # assert(unique2 == old_u2)
        
        logger.info("File: %s chromosome: %s Number of unique regions: %d" % (file1, chr,  len(unique1) ) )
        logger.info("File: %s chromosome: %s Number of unique regions: %d" % (file2, chr,  len(unique2) ) )
        logger.info("chromosome: %s Number of overlapping regions: %d" % (chr, len(overlap_list )) ) 
        
        write_unique_sites_per_chromosome_to_file(u1fid, chr, list1, unique1) 
        write_unique_sites_per_chromosome_to_file(u2fid, chr, list2, unique2) 
        
        for elem in overlap_list.keys(): 
            overlap_start =  overlap_list[elem][0]
            overlap_end =  overlap_list[elem][1]
            start_site1 = list1[elem[0] ][0][0]
            end_site1 = list1[elem[0] ][0][1]
            start_site2 = list2[elem[1] ][0][0]
            end_site2 = list2[elem[1] ][0][1]

            rest_line1 = list1[elem[0] ][1]
            rest_line2 = list2[elem[1] ][1]
            
            print >> cfid, chr, "\t", overlap_start,  overlap_end, "\t", start_site1, end_site1, "\t", start_site2, end_site2

            c1fid.write('%s\t%d\t%d\t%s\n' % (chr, start_site1, end_site1, rest_line1))


            c2fid.write('%s\t%d\t%d\t%s\n' % (chr, start_site2, end_site2, rest_line2))

            logger.debug('file: %s chromosome: %s common region: (%d, %d)' % (common1_file, chr, start_site1, end_site1))
            logger.debug('OVERLAPS WITH')
            logger.debug('file: %s chromosome: %s common region: (%d, %d)\n' % (common2_file, chr, start_site2, end_site2))
        logger.info('File: %s chromosome %s has %d common regions written to file: %s' % (file1, chr, len(overlap_list), common1_file))
        logger.info('File: %s chromosome %s has %d common regions written to file: %s' % (file2, chr, len(overlap_list), common2_file))
    u1fid.close()
    u2fid.close()
    cfid.close() 

    c1fid.close()
    c2fid.close() 

    fid1.close()
    fid2.close()
