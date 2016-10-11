#import os
import unittest
from unittest import TestCase
import numpy as np

#import knpackage.toolbox as kn
import gene_prioritization_toolbox as gptbx

class TestPerform_lasso_cv_regression(TestCase):

    def test_perform_lasso_cv_regression(self):
        spreadsheet = np.array([[0.8, 1.0, 0.4], [-0.3, 0.6, 0.1], [0.4, -0.1, 0.1], [0.4, -0.1, 0.1]])
        drug_response = np.array([[0.5, 0.7, 0.1]])
        expected_lasso = np.array([0.9990000000000001, 0.0, -0.0, -0.0])
        expected_predictor = np.array([0.49993333, 0.69973333, 0.10033333])
        lasso_array, lasso_predictor = gptbx.perform_lasso_cv_regression(spreadsheet, drug_response)

        for n in range(0, expected_lasso.size):
            self.assertAlmostEqual(expected_lasso[n], lasso_array[n], msg='lasso coefficent unexpected')

        for n in range(0, expected_predictor.size):
            self.assertAlmostEqual(expected_predictor[n], lasso_predictor[n], msg='lasso predictor unexpected')

if __name__ == '__main__':
    unittest.main()