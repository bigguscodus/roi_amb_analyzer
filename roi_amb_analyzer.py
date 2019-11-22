import argparse
from HLAparse import hla_parser
from HLAfilter_framer import hla_filter_framer
from HLAallele_solver import hla_allele_solver
from HLA_allele_writer import hla_allele_writer
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--imgt_file', default='hla.dat', type=str, help='Enter path to the database')
parser.add_argument('region', help="HLA prefix and gene. Example: HLA-A")
parser.add_argument('exon_numbers', help='Enter the numbers of exons, separated by space. Example: 2 3 4', nargs='+')
parser.add_argument('-a', '--amb',
                    help='Enter the path for ambiguous alleles output. Example ~/Documents/amb')
parser.add_argument('-u', '--unq',
                    help='Enter the path for nonambiguous alleles output. Example ~/Documents/uni')
args = parser.parse_args()
allele_exons_seqs = hla_parser(data=args.imgt_file, region=args.region)
data_frame = hla_filter_framer(allele_exon_dict=allele_exons_seqs,exon_numbers=args.exon_numbers)
allele_lists = hla_allele_solver(data_frame=data_frame)
hla_allele_writer(allele_lists=allele_lists,amb_allele_file=args.amb,uniq_allele_file=args.unq)

