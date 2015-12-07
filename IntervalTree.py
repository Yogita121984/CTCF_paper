# IntervalTree.py --- 
# 
# Description: Implementation of Interval Tree. 
# 
# Status: 
# Author: Vinay Jethava
# Created: Thu Jan 16 16:08:31 2014 (+0100)
# Last-Updated: 
#           By: 
#     Update #: 0
# 

# Change Log:
# 
# 

# Code:

import unittest
import sys

############################################################
## Need to have deeper recursion for tree construction
############################################################
MAX_RECURSION_LIMIT = 10000
sys.setrecursionlimit(MAX_RECURSION_LIMIT)

def median(x):
    return sorted(x)[len(x)//2]


class Interval:
    def __init__(self, l, r, id=None):
        self._left = l
        self._right = r
        self._id = id 

    def __str__(self):
        return ( '(' + str(self.left) + ', ' + str(self.right) + ')' )

    def __repr__(self): 
        return ( '(' + str(self.left) + ', ' + str(self.right) + ')' )
    
    @property
    def left(self):
        return self._left

    @property
    def right(self):
        return self._right

    @property
    def id(self):
        return self._id
    
    def find_overlap(self, other):
        b1 = self.left
        e1 = self.right
        b2 = other.left
        e2 = other.right
        overlap = None
        if ((e1 > b2) and (b1 < e2)) or ((b1 < e2) and (b2 < e1)): 
            overlap =  Interval(max(b1, b2), min(e1, e2))
        return overlap

    def has_overlap(self, other): 
        b1 = self.left
        e1 = self.right
        b2 = other.left
        e2 = other.right
        return ((e1 > b2) and (b1 < e2)) or ((b1 < e2) and (b2 < e1))
    
    def is_within(self, other): 
        return (self.left >= other.left) and (self.right <= other.right)

    def encloses(self, other):
        return other.is_within(self) 
    
    
        
        
class IntervalTree:
    def __init__(self, val, lc, rc, intervals): 
        self.key = val
        self.left_child = lc
        self.right_child = rc
        self.intervals = intervals
        # self._parent = None
        # self._sorted_b = sorted(self.intervals, key= lambda x : x.left)
        # self._sorted_e = sorted(self.intervals, key=lambda x: x.right, reverse=True)

    # ''' Sorted by start points in ascending order '''     
    # @property
    # def intervals_by_beginning(self):
    #     return self._sorted_b

    # ''' Sorted by end points in descending order '''
    # @property
    # def intervals_by_ending(self):
    #     return self._sorted_e 
     
    # @property
    # def parent(self):
    #     return self._parent

    # @parent.setter
    # def parent(self, par):
    #     self._parent = par 
        
    @classmethod
    def init_from_list(cls, interval_list, par=None): 
        if len(interval_list) == 0: 
            result = None
        else: 
            single_list = reduce(lambda x, y: x + [y.left , y.right ] , interval_list , [])
            median_val = median(single_list)
            left_list = filter(lambda x: x.right < median_val, interval_list)
            right_list = filter(lambda x: x.left > median_val, interval_list)
            middle_list = filter(lambda x: (x.right >= median_val) and (x.left <= median_val) , interval_list)
            result = cls(median_val, None, None , middle_list)
            if left_list is not None: 
                result.left_child = cls.init_from_list(left_list)
            if right_list is not None: 
                result.right_child = cls.init_from_list(right_list)
        return result 

    def overlapping_interval_search(self, Q):
        F = filter(lambda x: Q.has_overlap(x) , self.intervals)
        if (self.key >= Q.left) and (self.key <= Q.right):
            if self.left_child is not None:
                F = F + self.left_child.overlapping_interval_search(Q)
            if self.right_child is not None: 
                F = F  + self.right_child.overlapping_interval_search(Q)
        elif (self.key > Q.right) and (self.left_child is not None):
            F = F + self.left_child.overlapping_interval_search(Q)
        elif (self.key < Q.left) and (self.right_child is not None):
            F = F + self.right_child.overlapping_interval_search(Q)
        return F 

    def enclosing_interval_search(self, Q): 
        F = filter(lambda x: x.encloses(Q) , self.intervals)
        if (self.key >= Q.left) and (self.key <= Q.right):
            if self.left_child is not None:
                F = F + self.left_child.enclosing_interval_search(Q)
            if self.right_child is not None: 
                F = F  + self.right_child.enclosing_interval_search(Q)
        if (self.key > Q.right) and (self.left_child is not None):
            F = F + self.left_child.enclosing_interval_search(Q)
        elif (self.key < Q.left) and (self.right_child is not None):
            F = F + self.right_child.enclosing_interval_search(Q)
        return F 

class TestIntervalTree(unittest.TestCase): 
    def setUp(self):
        self.interval_list = [Interval(2, 3), Interval(0, 2),  Interval(0, 1), Interval(1.5, 2.5), Interval(4, 5), Interval(5, 6) ]
        self.interval_tree = IntervalTree.init_from_list(self.interval_list)
        
    def test_overlapping(self): 
        print '\ninterval_list' , self.interval_list
        Q1 = Interval(1.5 , 4.1)
        print 'Q', Q1
        print 'overlap', self.interval_tree.overlapping_interval_search(Q1)

    def test_enclosing(self): 
        print '\ninterval_list', self.interval_list
        Q2 = Interval(0.1, 0.2)
        print 'Q', Q2, 'enclosing', self.interval_tree.enclosing_interval_search(Q2)
        
if __name__=='__main__': 
    unittest.main()

