import argparse
from HLA_secondary_functions import hla_parser
from HLA_secondary_functions import hla_filter_framer
from HLA_secondary_functions import hla_allele_solver
from HLA_secondary_functions import hla_allele_writer

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--imgt_file', default='hla.dat', help='Enter path to the database')
parser.add_argument('-g', '--gene', required=True, help="HLA prefix and gene. Example: HLA-A")
parser.add_argument('-e', '--exon_numbers', required=True,
                    help='Enter the numbers of exons, separated by space. Example: 2 3 4',
                    nargs='+')
parser.add_argument('-a', '--amb',
                    help='Enter the path for ambiguous alleles output. Example ~/Documents/amb')
parser.add_argument('-u', '--unq',
                    help='Enter the path for nonambiguous alleles output. Example ~/Documents/uni')
parser.add_argument('-r', '--resolution', required=True,
                    help='Determines the level of detail of the alleles in a result according to nomenclature. '
                         'Levels are: 1,2,3,4')
args = parser.parse_args()
allele_exons_seqs = hla_parser(data=args.imgt_file, gene=args.gene)
data_frame = hla_filter_framer(allele_exon_dict=allele_exons_seqs, exon_numbers=args.exon_numbers)
allele_lists = hla_allele_solver(data_frame=data_frame, resolution=args.resolution)
hla_allele_writer(allele_lists=allele_lists, amb_allele_file=args.amb, uniq_allele_file=args.unq, gene=args.gene,
                  exon_numbers=args.exon_numbers, resolution=args.resolution)
