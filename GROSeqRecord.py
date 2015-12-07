
import re 

GROSEQ_DEFAULT_SEPARATOR="\t" 

############################################################
# Class for GROSeq records 
############################################################
class GROSeqRecord:
    '''
    GROSeqRecord - chromosome, peak_start, peak_end 
    '''
    def __init__(self, _chr, _ps, _pe, sep=GROSEQ_DEFAULT_SEPARATOR, strand=None): 
        self.separator = sep
        self.chromosome = _chr
        try:
            self.peak_start = int(_ps)
        except ValueError:
            self.peak_start = int(round(float(_ps)))
        try: 
            self.peak_end = int(_pe)
        except ValueError:
            self.peak_end = int(round(float(_ps))) 
        
        self.strand = strand
        

    ''' Class method that constructs GROSeqRecord from GROSeq string''' 
    @classmethod
    def read_record(cls, record, sep=GROSEQ_DEFAULT_SEPARATOR): 
        result = re.split('\s', record)
        if (len(result) == 6) and ((result[5] is '+') or  (result[5] is '-')):
            strand = result[5]
        else:
            strand = None
        obj = cls(result[0] , result[1], result[2], sep, strand)
        return obj 

    def __eq__(self, other):
        return ((self.chromosome, self.peak_start, self.peak_end) ==
                (other.chromosome, other.peak_start, other.peak_end))

    def __hash__(self):
        return hash((self.chromosome, self.peak_start, self.peak_end))

    def __cmp__(self, other): 
        first_cmp = cmp(self.chromosome, other.chromosome)
        if first_cmp is 0:
            return cmp(self.peak_start, other.peak_start)
        else:
            return first_cmp 
    def __repr__(self): 
        return self.separator.join([ self.chromosome, str(self.peak_start), str(self.peak_end) ])
    def __str__(self): 
        return self.separator.join([ self.chromosome, str(self.peak_start), str(self.peak_end) ])

############################################################
# Function from comparing two GROSeq string records
# @arg x: GROSeq record string 
# @arg y: GROSeq record string
# @return cmp_val: comparison between x and y 
#
# @note Same is achieved by GROSeqRecord class comparator
#       but sorting seems to be a bit faster using GROSeqRecord
# 
############################################################
def GROSeqStringCmp(x , y):
    rx = re.split('\s', x)
    ry = re.split('\s', y)
    first_cmp = cmp(rx[0], ry[0])
    second_cmp = cmp(int(rx[1]), int(ry[1]))
    if not (first_cmp == 0):
        return first_cmp
    else:
        return second_cmp 

