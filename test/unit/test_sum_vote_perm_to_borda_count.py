#import os
import unittest
from unittest import TestCase
import numpy as np

#import knpackage.toolbox as kn
import gene_prioritization_toolbox as gptbx

class TestSum_vote_perm_to_borda_count(TestCase):

    def test_sum_vote_perm_to_borda_count(self):

        pc_arr = np.array([0.162876, 0.800422, 0.080408, 0.842713, 0.750278, 0.969880])
        borda_count = np.zeros(pc_arr.size)
        expected_result = np.array([2, 4, 1, 5, 3, 6])
        borda_returned = gptbx.sum_vote_perm_to_borda_count(borda_count, pc_arr)
        for n in range(0, expected_result.size):
            self.assertEqual(expected_result[n], borda_returned[n], msg='borda count unexpected')

if __name__ == '__main__':
    unittest.main()