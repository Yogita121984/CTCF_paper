# test_gene_matrix_transposed.py --- 
# 
# Filename: test_gene_matrix_transposed.py
# Description: 
# Author: Vinay Jethava
# Created: Mon Nov  4 12:59:01 2013 (+0100)
# Version: 
# Last-Updated: Mon Nov  4 15:18:15 2013 (+0100)
#           By: Vinay Jethava
#     Update #: 32
# 

# Change Log:
# 
# 
# 
import unittest
from gene_matrix_transposed import GeneMatrixTransposed

class TestGeneMatrixTransposed(unittest.TestCase): 
    def setUp(self):
        self.test_file = 'example_gmt/example.gmt'
        self.gmt = GeneMatrixTransposed.read_gmt_file(self.test_file) 

    def testReadGMT(self): 
        self.out_gmt_file = 'example_gmt/copy_of_example.gmt'  
        self.gmt.write_gmt_file(self.out_gmt_file) 
        gmt2 = GeneMatrixTransposed.read_gmt_file(self.out_gmt_file)
        self.assertTrue(self.gmt == gmt2) 

    def testOverlap(self): 
        self.overlap_file = 'example_gmt/overlap.txt' 
        self.gmt.write_overlap_file(self.overlap_file) 

    def testMST(self): 
        a  = [ [0, 1, 4, 3 ], 
               [1, 0, 1, 2 ],
               [4, 1, 0, 5 ], 
               [3, 2, 5, 0 ] ]
        from UnionFind import UnionFind 
        from MinimumSpanningTree import MinimumSpanningTree 
        mst = MinimumSpanningTree(a) 
        self.assertTrue(mst == [(0, 1), (1, 2) , (1, 3)])
        
    def testAdjacencyList(self):
        threshold = 0.3
        self.gmt.write_adjacency_list('example_gmt/not_connected_0x3.txt', threshold, False)
        self.gmt.write_adjacency_list('example_gmt/connected_0x3.txt', threshold, True)
    
        # print 'Connected', len(self.gmt.get_adjacency_list(threshold, True))
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGeneMatrixTransposed)
    unittest.TextTestRunner(verbosity=3).run(suite)
