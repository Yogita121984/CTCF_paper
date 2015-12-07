#!/usr/bin/env python
#
# Description: This file counts the peaks occuring in a certain region on the chromosome.
# 
# Status: 
# Author: Vinay Jethava
# Created: Wed Jan 15 21:24:48 2014 (+0100)
# Last-Updated: 
#           By: 
#     Update #: 0
# 

# Change Log:
# 
# 

# Code:
import sys,logging,argparse
from GROSeqRecord import GROSeqRecord
from update_progress import update_progress
from extract_common import read_file_to_map
from logger import logger,set_verbosity  
from IntervalTree import Interval, IntervalTree


############################################################
## Default options - can be changed at command line 
############################################################

## strand specific. 
STRAND_SPECIFIC = True

## Only write regions with non-zero number of peaks
SKIP_ZERO_COUNTS = True

## How much logging do you want? 0 = minimal, 2 = verbose
VERBOSITY = 1 

############################################################
# Increase recursion level
############################################################
MAX_RECURSION_LEVEL = 20000
sys.setrecursionlimit(MAX_RECURSION_LEVEL) 



############################################################
# Count peaks in strand_specific manner
############################################################
def count_peaks_in_regions(chromosome_region_map, groseq_peaks):
    ''' Count peaks falling in different chromosome regions '''
    logger.info('Reading groseq data')
    max_iter = len(groseq_peaks)
    curr_iter = 0
    result = dict()
    chr_interval_trees = dict()
    ############################################################
    ## Initialize result
    ############################################################
    for chr in chromosome_region_map.keys():
        num_regions = len(chromosome_region_map[chr])
        result[chr] = ([0]*num_regions, [0]*num_regions)
        ############################################################
        ## Improved timing using interval trees.
        ############################################################
        interval_list = []
        idx = 0        
        logger.info('chromosome ' + chr + ' has regions: ' + str(num_regions))
        for region in chromosome_region_map[chr]:
            update_progress(idx, 'Reading regions' , num_regions, 1000)
            r_start = region[0][0]
            r_end = region[0][1]
            interval_list = interval_list + [ Interval(r_start, r_end, idx) ]
            idx = idx + 1
        logger.info('Creating IntervalTree for chromosome' + chr) 
        chr_interval_trees[chr] = IntervalTree.init_from_list(interval_list)
        
    ############################################################
    ## Processing GROSeq begins here
    ############################################################
    for peak in groseq_peaks:
        peak = peak.strip()
        ############################################################
        ## progress line below - can be removed.
        ############################################################
        curr_iter = curr_iter + 1
        update_progress(curr_iter, ('Reading GROSeq peaks file'), max_iter)

        curr_record = GROSeqRecord.read_record(peak)
        curr_chromosome  = curr_record.chromosome
        if curr_chromosome in chromosome_region_map:
            regions_on_curr_chromosome = chromosome_region_map[curr_chromosome]
            peak_interval = Interval(curr_record.peak_start, curr_record.peak_end )
            enclosing_intervals = chr_interval_trees[curr_chromosome].enclosing_interval_search(peak_interval)
            ############################################################
            ## update the counts
            ############################################################
            curr_count = result[curr_chromosome]
            curr_plus= curr_count[0]
            curr_minus = curr_count[1]
            for interval in enclosing_intervals:
                logger.debug('peak: ' + str(peak_interval) + ' on ' + curr_chromosome + ' ' + str(interval) + ' strand: ' +  curr_record.strand)
                idx = interval.id
                if curr_record.strand is '+':
                    curr_plus[idx] = curr_plus[idx] + 1
                elif curr_record.strand is '-':
                    curr_minus[idx] = curr_minus[idx] + 1
            result[curr_chromosome] = (curr_plus, curr_minus)
    return result

############################################################
# This function does all the lifting
############################################################
def main(chromosome_region_file, groseq_peak_file, output_file, strand_specific=STRAND_SPECIFIC, skip_zero_counts=SKIP_ZERO_COUNTS, separator='\t'): 
    ''' Main function that operates on the files 

    @arg chromosome_region_file BED file containing list of chromosome regions 
    @arg groseq_peak_file       GROSEQ file containing groseq peaks
    @arg output_file            output file (output_file+ and output_file- if strand_specific is True)
    
    '''
    logger.info('bed_file: ' + chromosome_region_file)
    logger.info('groseq_peak_file: ' + groseq_peak_file)
    logger.info('output_file: ' + output_file) 
    logger.info('strand_specific: ' + str(strand_specific) ) 
    logger.info('skip_zero_counts: ' +  str(skip_zero_counts) )
    with open(chromosome_region_file , 'r') as chr_fid:
        ############################################################
        ## Construct the chromosome chr -> [(region_start, region_end)*]
        ############################################################
        logger.info('Reading chromosome region file') 
        chromosome_region_map = read_file_to_map(chr_fid, separator)
        with open(groseq_peak_file, 'r') as groseq_fid:
            ## Read next peak from groseq file
            groseq_lines = groseq_fid.readlines()
            result = count_peaks_in_regions(chromosome_region_map, groseq_lines)                

            if strand_specific is True:
                plus_fid = open((output_file + '+'), 'w')
                minus_fid = open((output_file + '-'), 'w')
                for chr in chromosome_region_map.keys():
                    curr_regions = chromosome_region_map[chr]
                    curr_num_peaks = result[chr]
                    n = len(curr_regions)
                    r_count_p = 0
                    r_count_n = 0                         
                    logger.info('Writing strand specific output for ' + chr + ' having ' + str(n) + ' regions')
                    for j in range(n):    
                        update_progress(j, 'Writing output for chromosome ' + chr, n)
                        region_start = curr_regions[j][0][0]
                        region_end = curr_regions[j][0][1]
                        num_peaks_plus = curr_num_peaks[0][j]
                        num_peaks_minus = curr_num_peaks[1][j]

                        if (not skip_zero_counts) or (num_peaks_plus > 0):
                            r_count_p = r_count_p + 1 
                            print >> plus_fid, separator.join([chr, str(region_start), str(region_end), str(num_peaks_plus)])
                        if (not skip_zero_counts) or (num_peaks_minus > 0):
                            r_count_n = r_count_n + 1 
                            print >> minus_fid, separator.join([chr, str(region_start), str(region_end), str(num_peaks_minus)])
                    logger.info('Wrote ' + str(r_count_p) + ' of '+ str(n) + ' regions for ' + chr + ' with non-zero peaks on strand: + ' )
                    logger.info('Wrote ' + str(r_count_n) + ' of '+ str(n) + ' regions for ' + chr + ' with non-zero peaks on strand: - ' )
                plus_fid.close()
                minus_fid.close()
            else:
                with open(output_file, 'w') as o_fid:
                    for chr in chromosome_region_map.keys():
                        curr_regions = chromosome_region_map[chr]
                        n = len(curr_regions)
                        regions_written = 0 
                        for j in range(n):
                            update_progress(j, 'Writing output for chromosome ' + chr, n)
                            region_start = curr_regions[j][0][0]
                            region_end = curr_regions[j][0][1]
                            num_peaks = result[chr][0][j] + result[chr][1][j]
                            if (not skip_zero_counts) or (num_peaks > 0): 
                                regions_written = regions_written + 1
                                print >>o_fid, separator.join([chr, str(region_start), str(region_end), str(num_peaks)])
                        logger.info('Finished writing ' + str(regions_written) + ' of '+ str(n) + ' regions for chromosome ' + chr)
                        

############################################################
# TODO: Argument parser
############################################################

# parser.add_argument('-p', metavar='UNIQUE_PREFIX', type=str, help='Prefix appended to unique reads in file (default: u_)', default='u_') 
# parser.add_argument('-o', metavar='OVERLAP_FILE', type=str, help='File containing region overlaps (default: overlap.txt)', default='overlap.txt')



if __name__=='__main__':
    ############################################################
    ## Default settings
    ############################################################
    skip_zero_counts = SKIP_ZERO_COUNTS
    strand_specific = STRAND_SPECIFIC
    logger.setLevel(logging.INFO) 
    
    parser = argparse.ArgumentParser(description='This script counts the peaks occuring in a certain region on the chromosome.')
    parser.add_argument('bed_file', metavar='CTCF_BED_FILE', type=str, help='CTCF file containing chromosome regions')
    parser.add_argument('groseq_file', metavar='GROSEQ_PEAK_FILE', type=str, help='GROSEQ file containing peak information') 
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str, help='In strand specific mode, the program writes OUTPUT_FILE+ and OUTPUT_FILE- (for + and - strand resp.) else (if merge_strands option is given) generates single OUTPUT_FILE with counts for the two strands added together')  
    parser.add_argument('--keep_zero_counts',  action='store_true', help='Print regions with zero peaks')
    parser.add_argument('--merge_strands' ,  action='store_true', help='Merge counts from the two strands (+ and -)') 
    parser.add_argument('-v', '--verbosity', metavar='VERBOSITY', help="Increase verbosity level [0, 1, 2] (default: 1)", type=int, default=VERBOSITY) 

    args = parser.parse_args(sys.argv[1:])

    bed_file = args.bed_file
    groseq_file = args.groseq_file
    output_file = args.output_file
    
    if args.keep_zero_counts:
        skip_zero_counts = False
    if args.merge_strands:
        strand_specific = False

    set_verbosity(args.verbosity) 
    # if args.verbosity == 0:
    #     logger.setLevel(logging.ERROR)
    # elif args.verbosity == 2:
    #     logger.setLevel(logging.DEBUG) 
    # else:
    #     logger.setLevel(logging.INFO)
    
    ############################################################
    # calling main function
    ############################################################
    main(bed_file, groseq_file, output_file, strand_specific, skip_zero_counts)


			

		

		
