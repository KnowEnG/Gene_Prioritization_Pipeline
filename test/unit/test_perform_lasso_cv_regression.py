#import os
import unittest
from unittest import TestCase
import numpy as np

#import knpackage.toolbox as kn
import gene_prioritization_toolbox as gptbx

class TestPerform_lasso_cv_regression(TestCase):

    def test_perform_lasso_cv_regression(self):
        samp_feat_matrix = np.array([[0.8, 1.0, 0.4], [-0.3, 0.6, 0.1], [0.4, -0.1, 0.1], [0.4, -0.1, 0.1]])
        response_array = np.array([[0.5, 0.7, 0.1]])
        expected_lasso = np.array([0.9990000000000001, 0.0, -0.0, -0.0])
        expected_predictor = np.array([0.49993333, 0.69973333, 0.10033333])
        lasso_array = gptbx.perform_lasso_cv_regression(samp_feat_matrix, response_array)

        for n in range(0, expected_lasso.size):
            self.assertAlmostEqual(expected_lasso[n], lasso_array[n], msg='lasso coefficent unexpected')

if __name__ == '__main__':
    unittest.main()