from functions import run
import argparse

parser = argparse.ArgumentParser(description='Parser generating kmers from pdb files of RoseTTAFold2NA')
parser.add_argument('--ncores', type=int, default=1, help='The number of cores used for extracting kmers from pdb files')
parser.add_argument('--threshold', type=float, default=0.45, help='Threshold for distance between nucleic acid and amino acids for counting contact')
parser.add_argument('--scheme', type=str, choices=['ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'],
                    default='ca', help="Type of counting distance between nucleic acid and amino acid. "
                                       "Allowed values: 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'. "
                                       "For more information, visit: https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_contacts.html")
parser.add_argument('--dataset', type=str, choices=['sample', 'full'], default='sample', help='Sample or full dataset')

args = parser.parse_args()

ncores = args.ncores
threshold = args.threshold
scheme = args.scheme
dataset = args.dataset

if __name__ == '__main__':

    run(ncores, threshold, scheme, dataset)