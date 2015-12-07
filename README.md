# Introduction 

This repository contains generic scripts used for the analysis of ChipSeq and GROseq data used in CTCF paper.


* extract_common.py 
* gmt_main.py
* count_peaks.py 


## extract_common.py

This script extracts the common reads from two .bed files. This script is used to categorize the genomic distribution of CTCF binding sites.


     Usage: extract_common.py [-h] [-c COMMON_FILE] [-p UNIQUE_PREFIX] FILE1 FILE2

#### Input file format: 
   
    <chromosome> <start_site> <end_site> <REST_OF_LINE IN ORIGINAL FILE> 

#### Output files: 

- u_FILE1: Unique reads in FILE1 (full line copied from FILE1)
- u_FILE2: Unique reads in FILE2 (full line copied from FILE1)
- co_FILE1 and co_FILE2: Overlapping reads (same lines overlap) between FILE1 and FILE2
- common.txt: Line n shows positions at which co_FILE1 line n and co_FILE2 line n overlap. 

#### Format of common.txt: 
 
    <chromosome> <overlapping start> <overlapping end> <FILE1 region start> <FILE1 region end> <FILE2 region start> <FILE2 region end> 

#### Example output. 

For example, run the following commands: 

    $ cd example_extract_common/
    $ python ../extract_common.py TH19.bed TH20.bed

## gmt_main.py 

This code analyses GMT (GeneMatrixTransposed) files. 
     
     Usage: gmt_main.py -i INPUT_FILE -t THRESHOLD 

See gene_matrix_transposed.py for the class description. 

#### Output files: 

- overlap.csv
- gsID.csv 
- adj_list.csv 

See example_gmt/ for examples for generated files. 

#### Example run

    $ python gmt_main.py -i example_gmt/example.gmt -t 0.2 

## count_peaks.py

This file counts number peaks at a chromosome location from GRO-Seq
.bed file.

	Usage: $ count_peaks.py INPUT_FILE OUTPUT_FILE

The input file format is:

	CHROMOSOME PEAK_START PEAK_END n X STRAND(+/-)

and the corresponding output file format is:

	CHROMOSOME PEAK_START PEAK_END #TIMES_IN_FILE STRAND

In case of conflict in strands or missing strand, the strand
information is not given. 
	
#### Example 
Input file: 

	chr1	10554	10555	n	2	-
	chr1	10554	10555	n	2	-
	chr4	13365	13366	n	0	-
	chr1	13365	13366	n	0	-
	chr1	13365	13366	n	0	-
	chr3	14022	14023	n	0	-
	chr2	14022	14023	n	0	-
	chr1	14022	14023	n	0	-
	chr1	14022	14023	n	0	-
	chr1	14022	14023	n	0	-
	chr1	14022	14023	n	0	-
	chr1	14022	14023	n	0	-

The corresponding output is:

	chr1	10554	10555	2	-
	chr1	13365	13366	2	-
	chr1	14022	14023	5	-
	chr2	14022	14023	1	-
	chr3	14022	14023	1	-
	chr4	13365	13366	1	-

For testing, use the command

	python count_peaks.py example_peak_count/input.groseq example_peak_count/test.out

## count\_peaks\_in\_regions.py

This script takes a BED file (chromosome regions) and a GROSEQ peaks
file and counts the number of peaks (per strand) falling in each
chromosome region.

    Usage: count_peaks_in_region.py BED_FILE GROSEQ_FILE OUTPUT_FILE

Creates OUTPUT\_FILE+ and OUTPUT\_FILE- corresponding to two strands.

### Global options
    ### Setting this option to False counts the peaks irrespective of strand
    strand_specific = True

    ### Commenting this line would output logging messages
    logger.setLevel(logging.ERROR)

#### Example 

    ./count_peaks_in_region.py example_count_peaks_in_region/input.regions example_count_peaks_in_region/input2.groseq example_count_peaks_in_region/test.out
    
The above command generates `test.out+` and `test.out-` files in  `example_count_peaks_in_region` directory. 
