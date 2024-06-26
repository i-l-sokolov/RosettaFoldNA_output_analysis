from functions import *
from scipy.stats import chi2_contingency
import pandas as pd
from pandas_parallel_apply import DataFrameParallel, SeriesParallel
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
import warnings

parser = argparse.ArgumentParser(description='Parser generating kmers from pdb files of RoseTTAFold2NA')
parser.add_argument('--ncores', type=int, default=1, help='The number of cores used for exctracting kmers from pdb files')
parser.add_argument('--seed', type=int, default=42, help='Seed for generating random kmers')
parser.add_argument('--threshold', type=float, default=0.45, help='Threshold for distance between nucleic acid and amino acids for counting contact')
parser.add_argument('--scheme', type=str, choices=['ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'],
                    default='ca', help="Type of counting distance between nucleic acid and amino acid. "
                                       "Allowed values: 'ca', 'closest', 'closest-heavy', 'sidechain', 'sidechain-heavy'. "
                                       "For more information, visit: https://mdtraj.org/1.9.4/api/generated/mdtraj.compute_contacts.html")

args = parser.parse_args()

ncores = args.ncores
seed = args.seed
threshold = args.threshold
scheme = args.scheme

if __name__ == '__main__':
    #Reading the input dataframe with information about pdb files
    df_input = pd.read_csv('../data/df_input.csv')

    # Exctracting kmers from pdb files
    print('Extracting kmers with {} cores...'.format(ncores))
    if ncores == 1:
        df_input['extracted_kmers'] = df_input['pdb_file'].apply(lambda x : get_kmer(x, threshold=threshold, scheme=scheme))
    else:
        #Creating function as SeriesParallel doesn't take lambda function
        def get_kmer_arguments(x):
            return get_kmer(x, threshold=threshold, scheme=scheme)
        #Shutting down warnings as SeriesParallel, probably, uses np.split for split pandas DataFrame
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            df_input['extracted_kmers'] = SeriesParallel(df_input['pdb_file'], n_cores=ncores, pbar=False).apply(get_kmer_arguments)
    print('Kmers have been extracted')

    #Generation random set of kmers
    df_input['rand_kmer'] = df_input['DNA'].apply(lambda x : random_choice(x, seed=seed))
    print('Random kmers have been generated')

    #Calculating score for exctracted and random generated kmers
    df_input['RF_gen_score'] = df_input.apply(lambda x : score(x, 'extracted_kmers'), axis=1)
    df_input['random_gen_score'] = df_input.apply(lambda x : score(x, 'rand_kmer'), axis=1)
    print('Scores for extracted and random kmers have been calculated')

    #Creating order of modes for resulting picture
    order = sorted([x for x in df_input['mode'].unique() if x[-2:] not in ['_r','_1']]) + \
    sorted([x for x in df_input['mode'].unique() if x[-2:] in ['_r','_1']])

    #Saving DataFrame with results
    df_input.to_csv('../results/df_res.csv', index=False)
    print('Dataframe with results has been saved to results/df_res.csv')

    #Creating resulting picture
    g3 = sns.FacetGrid(pd.melt(df_input[['random_gen_score','RF_gen_score','mode','reverse']].rename(columns={'RF_gen_score' : 'Rosetta_sc', 'random_gen_score' : 'rand_sc'}), id_vars=['mode','reverse'], var_name='score', value_name='results'), col='reverse', row='score')
    g3.map_dataframe(sns.countplot, x='results', hue='mode', order=[True,False],
                    hue_order=order,
                    palette = 'tab10')
    plt.legend(bbox_to_anchor=(1.6,1),loc='center right')
    g3.savefig('../results/rosetta.png', dpi=1000)
    g3.savefig('../results/rosetta.pdf', format='pdf', dpi=1000)
    g3.savefig('../results/rosetta.svg', format='svg', dpi=1000)
    print('Figure has been saved to results/rosetta in png, pdf and svg formats')

    #Calculating Chi Square p values and saving them as new DataFrame
    df_input[['random_gen_score','RF_gen_score','mode']].groupby('mode').apply(
        lambda x : chi2_contingency(pd.crosstab(x['random_gen_score'],x['RF_gen_score']))[1], include_groups=False
    ).reset_index().rename(columns={0:'Chi Square Test p-values', 'mode' : 'Mode'}).to_csv('../results/chi2_contingency.csv', index=False)
    print('Chi Square Test p values have been saved to results/chi2_contingency.csv')