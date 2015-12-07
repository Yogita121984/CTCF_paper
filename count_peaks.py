#!/usr/bin/env python 
# count_peaks.py --- 
# 
# USAGE: python count_peaks.py  INPUT_FILE OUTPUT_FILE
#  

# Filename: count_peaks.py
# Description: This file counts the number of peaks in GRO-seq .bed file. 
# Author: Vinay Jethava

# Input format (See example_peak_count.groseq): 
# 
# chr1	10554	10555	n	2	-
# chr1	10554	10555	n	2	-
# chr1	13365	13366	n	0	-
# chr1	13365	13366	n	0	-
# chr1	13365	13366	n	0	-


# Output format (See peak_count.output): 
#
# chr1 10554 10555 2 - 
# chr1 13365 13366 3 - 
# 
 
# Change Log:
# 14-Jan-2014   Vinay Jethava
#    Updated to using the comparator GROSeqStringCmp() instead of the class
# 5-Jan-2014    Vinay Jethava  
#    Updated version - strand conflict detection. Sorted output. 
# 5-Jan-2014    Vinay Jethava  
#    initial version with input/output specification
#  

import sys
import re 
import math

from logger import logger
from GROSeqRecord import GROSeqRecord,GROSeqStringCmp
from update_progress import update_progress, PROGRESS_DELAY 

############################################################
# Global parameters
############################################################
''' Separator for output file formatting '''
SEPARATOR = '\t'


############################################################
# Peak counter
############################################################    
def count_peaks(input_file, output_file): 
    '''
    Top-level function which counts peaks from GROSeq files

    @param input_file: Input GROSeq .bed file
    @param output_file: File to which write the output
    '''
    with open(input_file, 'r') as fid: 
        line_count = 0
        lines = fid.readlines() 
        res_dict = dict()
        total_lines = len(lines)   
        unreliable_strands_flag = False
        ############################################################
        # Read input file
        ############################################################
        for line_orig in lines:
            line_count = line_count + 1
            ############################################################
            #  Updating progress counter - can be safely removed
            ############################################################
            if line_count % PROGRESS_DELAY == 0:
                progress = line_count * 1.0 / total_lines 
                update_progress(progress, ('Reading ' + input_file)) 
            ############################################################
            # Main function body
            ############################################################    
            line_orig = line_orig.strip()
            parts = re.split('\s', line_orig)
            ### do not put strand unless specified
            strand = ''             
            if len(parts) < 3: 
                print >> sys.stderr,  line_count,  line_orig, 'Incomplete line - skipping'
                continue
            elif len(parts) < 6: 
                unreliable_strands_flag = True
                print >> sys.stderr, line_count, line_orig, 'missing strand - ignoring strands in output!'              
            else: 
                strand = parts[5]
            ### constructing basic line - faster using str as key!
            line = SEPARATOR.join(parts[0:3])
            if line in res_dict:
                prev_count = res_dict[line][0]
                prev_strand = res_dict[line][1]
                ############################################################
                ## THIS LINE NEEDS TO BE REMOVED IF WISH TO SEE CONFLICTS
                ############################################################
                if (not unreliable_strands_flag) and (strand is not prev_strand):
                    print >> sys.stderr, "\n", line_count, line, ' strand: ', strand, ' prev_strand: ' , prev_strand
                    unreliable_strands_flag = True
                res_dict[line] = (prev_count + 1, prev_strand)
            else:
                res_dict[line] = (1, strand)
        ############################################################
        # Write output file
        ############################################################
        output_line_count = len(res_dict)
        curr_oline = 0
        with open(output_file, 'w') as ofid:
            ############################################################
            ### sorting based on GROSeqRecord structure. 
            ############################################################
            print >> sys.stderr, "\nSorting", output_line_count, "records for printing"             
            records = [ GROSeqRecord.read_record(x) for x in res_dict.iterkeys() ]
            for record in sorted(records):                
                key = str(record)
                curr_oline = curr_oline + 1 
                ############################################################
                # Progress counter - can be safely removed
                ############################################################
                if curr_oline % PROGRESS_DELAY == 0: 
                    progress = curr_oline * 1.0 / output_line_count
                    update_progress(progress, ('Writing ' + output_file))  
                ############################################################
                # Writing output line 
                ############################################################    
                oline = str(key) + SEPARATOR + str(res_dict[key][0])
                if unreliable_strands_flag is False:
                    oline = oline + SEPARATOR  + res_dict[key][1]
                print >> ofid , oline


if __name__=='__main__':
    if len(sys.argv) < 2: 
        sys.exit('Usage: count_peaks.py INPUT_FILE OUTPUT_FILE')
    input_file = sys.argv[1]
    output_file = sys.argv[2]     
    count_peaks(input_file, output_file) 
     

        
