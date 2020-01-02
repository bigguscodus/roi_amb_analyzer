import unittest
import pandas as pd
import HLA_secondary_functions as hla

test_dict = {'HLA-TEST1*01:01': {1: 'AAAA', 2: 'TTTT'},
             'HLA-TEST1*01:02': {1: 'GGGG', 2: 'CCCC'},
             'HLA-TEST1*01:03': {1: 'AAAA', 2: 'TTTT'},
             'HLA-TEST1*01:04': {1: 'GGG', 2: 'TTT'},
             'HLA-TEST1*02:05': {1: 'AAAA'},
             'HLA-TEST1*02:01:01': {1: 'GGG', 2: 'CCC', 3: 'AAAA'},
             'HLA-TEST1*02:01:02': {1: 'GGG', 2: 'CCC', 3: 'AAAA'},
             }
test_df = pd.DataFrame(test_dict)
test_df = test_df.loc[[1, 2, 3], ].reset_index().drop(labels='index', axis=1)
test_tuple = ([{'HLA-TEST1*01:01', 'HLA-TEST1*01:03', 'HLA-TEST1*02:05'}, {'HLA-TEST1*02:01:01', 'HLA-TEST1*02:01:02'}],
             ['HLA-TEST1*01:02', 'HLA-TEST1*01:04'])

class MyTestCase(unittest.TestCase):
    def test_output_types(self):
        self.assertIsInstance(hla.hla_parser(data = '/home/mikka/PycharmProjects/parseq_lab/hla.dat', gene = 'HLA-A'), dict)
        self.assertIsInstance(hla.hla_filter_framer(test_dict, exon_numbers=['1', '2', '3']), pd.DataFrame)
        self.assertIsInstance(hla.hla_allele_solver(test_df, resolution='4'), tuple)
        self.assertIsNone(hla.hla_allele_writer(test_tuple,amb_allele_file='del',uniq_allele_file=None,gene='HLA-A',exon_numbers=[2,4],resolution='4'))

    def test_output_values(self):
        self.assertEqual(
            ([{'HLA-TEST1*01:01', 'HLA-TEST1*01:03', 'HLA-TEST1*02:05'}, {'HLA-TEST1*02:01:01', 'HLA-TEST1*02:01:02'}],
             ['HLA-TEST1*01:02', 'HLA-TEST1*01:04']),
            hla.hla_allele_solver(test_df, resolution='4'))
        self.assertEqual(
            ([{'HLA-TEST1*01:01', 'HLA-TEST1*01:03', 'HLA-TEST1*02:05'}, {'HLA-TEST1*02:01:01', 'HLA-TEST1*02:01:02'}],
             ['HLA-TEST1*01:02', 'HLA-TEST1*01:04']),
            hla.hla_allele_solver(test_df, resolution='3'))
        self.assertEqual(
            ([{'HLA-TEST1*01:01', 'HLA-TEST1*01:03', 'HLA-TEST1*02:05'}],
             ['HLA-TEST1*01:02', 'HLA-TEST1*01:04']),
            hla.hla_allele_solver(test_df, resolution='2'))
        self.assertEqual(
            ([{'HLA-TEST1*01', 'HLA-TEST1*02'}],
             ['HLA-TEST1*01:02', 'HLA-TEST1*01:04']),
            hla.hla_allele_solver(test_df, resolution='1'))
    def test_input_errors(self):
        self.assertIsNone(hla.hla_parser(gene='HLA-X'))
        self.assertRaises(FileNotFoundError, hla.hla_parser(data='hla_1', gene='HLA-A'))
if __name__ == '__main__':
    unittest.main()
