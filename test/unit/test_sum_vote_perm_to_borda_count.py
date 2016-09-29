#import os
import unittest
from unittest import TestCase
import numpy as np

#import knpackage.toolbox as kn
import gene_prioritization_toolbox as gptbx

class TestSum_vote_perm_to_borda_count(TestCase):

    def test_sum_vote_perm_to_borda_count(self):
        borda_count = np.array([0,0,0,0,0])
        vote_rank = np.array([0,3,1,2])
        vote_perm = np.array([2,1,4,3])
        borda_expected = np.array([0,1,4,2,3])
        borda_returned = gptbx.sum_vote_perm_to_borda_count(borda_count, vote_rank, vote_perm)
        for n in range(0, borda_expected.size):
            self.assertEqual(borda_expected[n], borda_returned[n], msg='borda count unexpected')

        vote_rank = np.array([1,0,2,3])
        vote_perm = np.array([1,2,3,0])
        borda_expected = np.array([1,4,8,4,3])
        borda_returned = gptbx.sum_vote_perm_to_borda_count(borda_count, vote_rank, vote_perm)
        for n in range(0, borda_expected.size):
            self.assertEqual(borda_expected[n], borda_returned[n], msg='borda count unexpected')

        vote_rank = np.array([2,0,1,3])
        vote_perm = np.array([0,3,1,2])
        borda_expected = np.array([3,7,9,8,3])
        borda_returned = gptbx.sum_vote_perm_to_borda_count(borda_count, vote_rank, vote_perm)
        for n in range(0, borda_expected.size):
            self.assertEqual(borda_expected[n], borda_returned[n], msg='borda count unexpected')

if __name__ == '__main__':
    unittest.main()