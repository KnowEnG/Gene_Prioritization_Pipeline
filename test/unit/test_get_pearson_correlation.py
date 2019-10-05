import unittest
from unittest import TestCase
import numpy as np

import kngeneprioritization.gene_prioritization_toolbox as gpt

class TestGet_pearson_correlation(TestCase):
    """  predicted correlation vs gene_prioritization_toolbox.get_correlation for 'correlation_measure' = 'pearson'  """
    def test_get_pearson_correlation(self):
        """ (Note that this is a test of expected typical input and does not test extreme inputs or all possible inputs)
        1) data with same slope will have positive correlation and equal 1.0
        2) data with opposite slope will have negative correlation and equal -1.0
        3) data with orthogonal slope will equal zero
        4) there will be no non-numeric values returned from expected numeric input
        """
        run_parameters = {'correlation_measure': 'pearson'}
        n_test_rows = 20
        n_test_cols = 11
        positive_correlation_member = 3
        positive_correlation_value = 1.0
        negative_correlation_member = 4
        negative_correlation_value = -1.0
        no_correlation_member = 5
        no_correlation_value = 0.0

        spreadsheet_mat = np.random.random((n_test_rows, n_test_cols))
        drug_response = spreadsheet_mat[positive_correlation_member,:]

        spreadsheet_mat[negative_correlation_member,:] = drug_response * (-1) + 1
        spreadsheet_mat[no_correlation_member] =\
                spreadsheet_mat[negative_correlation_member,:] + spreadsheet_mat[positive_correlation_member,:]

        corr_arr = gpt.get_correlation(spreadsheet_mat, drug_response, run_parameters)

        self.assertAlmostEqual(corr_arr[positive_correlation_member], positive_correlation_value,
                               msg='Pearson correlation error')
        self.assertAlmostEqual(corr_arr[negative_correlation_member], negative_correlation_value,
                               msg='Pearson correlation error')
        self.assertAlmostEqual(corr_arr[no_correlation_member], no_correlation_value,
                               msg='Pearson correlation error')

        for correlation_n in corr_arr:
            self.assertTrue(np.isfinite(correlation_n), msg='Pearson Correlation Array is Not Finite')

if __name__ == '__main__':
    unittest.main()
