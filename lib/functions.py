import numpy as np
from biopandas.pdb import PandasPdb
import mdtraj as md
import itertools
import random


def random_choice(dna, seed):
    """
    Generates a random k-mer from the given DNA sequence.
    The function replaces a random number of elements in the DNA sequence with the symbol "_".

    :param dna: str, the DNA sequence from which to generate the k-mer
    :param seed: int, the seed for the random number generator to ensure reproducibility
    :return: str, the DNA sequence with random elements replaced by "_"
    """
    random.seed(seed)
    dna = list(dna)
    k = random.choice(range(len(dna)))
    reps = random.sample(range(0,len(dna)), k=k)
    for i in reps:
        dna[i] = "_"
    return "".join(dna)


def score(df, column):
    """
    Evaluates the number of underscore ("_") symbols in specific segments of a k-mer column
    from a given DataFrame, based on the mode of the DataFrame.

    :param df: pandas.DataFrame, the DataFrame containing the k-mer sequences and additional data.
    :param column: str, the column name in the DataFrame containing the k-mer sequence to be evaluated.
    :return: bool, True if the count of underscores in the first segment is greater or lesser
                   than the second segment based on the mode, otherwise False.
    """
    # Calculate the length for k-mer segments
    ln = len(df['kmer1']) + len(df['kmer2']) + int(df['dist'].split('N')[0])

    # Get the k-mer sequence from the specified column, removing any spaces
    kmers = df[column].replace(" ","")

    # Check the mode and perform the comparison based on the suffix of the mode. Reverse case should end with either _r or _1. The length of flanks is 5. Comparison only takes place in the region of DNA that potentially binds to TFs, without considering flanks or the distance between binding sites
    if df['mode'][-2:] in ['_r','_1']:
        return kmers[5:ln+5].count("_") > kmers[-ln-5:-5].count("_")
    else:
        return kmers[5:ln+5].count("_") < kmers[-ln-5:-5].count("_")


def get_mask(pdb, threshold, scheme):
    """
    Computes the contact mask for a given protein-DNA complex from a PDB file.

    Parameters:
    pdb (str): The filename of the PDB file to be read.
    threshold (float): The distance threshold to determine if a contact exists.
    scheme (str): The scheme used for computing contacts. This could be any valid scheme recognized by the md.compute_contacts function from the MDTraj library.

    Returns:
    numpy.ndarray: A boolean array indicating whether the minimum contact distance for each DNA residue is below the specified threshold.
    """

    # Reading pdb file
    pdb_file = PandasPdb()
    pdb_file.read_pdb('../data/commit_6a1b3d_pdbs/' + pdb)

    # Creating dict with indexes of chains
    tdict = dict(
        pdb_file.amino3to1(fillna="").reset_index(drop=True).groupby(by='chain_id').apply(lambda x: x.index.tolist(), include_groups=False))

    # getting pairs
    pairs_A_C = list(itertools.product(tdict['A'], tdict['C']))
    pairs_B_C = list(itertools.product(tdict['B'], tdict['C']))
    pairs_A_D = list(itertools.product(tdict['A'], tdict['D']))
    pairs_B_D = list(itertools.product(tdict['B'], tdict['D']))

    file = md.load('../data/commit_6a1b3d_pdbs/' + pdb)

    # Computing contacts
    contacts_A_C = md.compute_contacts(file, pairs_A_C, scheme=scheme)
    contacts_B_C = md.compute_contacts(file, pairs_B_C, scheme=scheme)
    contacts_A_D = md.compute_contacts(file, pairs_A_D, scheme=scheme)
    contacts_B_D = md.compute_contacts(file, pairs_B_D, scheme=scheme)

    # Creating arrays from contacts
    arr_A_C = contacts_A_C[0][0].reshape(len(tdict['A']), len(tdict['C'])).T
    arr_B_C = contacts_B_C[0][0].reshape(len(tdict['B']), len(tdict['C'])).T
    arr_A_D = contacts_A_D[0][0].reshape(len(tdict['A']), len(tdict['D'])).T
    arr_B_D = contacts_B_D[0][0].reshape(len(tdict['B']), len(tdict['D'])).T

    # Creating arrays for every DNA strand. For chain D the flip is applied because D is reverse compliment of the C
    df_C = np.concatenate([arr_A_C, arr_B_C], axis=1)
    df_D = np.flip(np.concatenate([arr_A_D, arr_B_D], axis=1), 0)

    # Stacking two arrays and choosing minimal through the first axe
    df = np.stack([df_C, df_D]).min(axis=0)

    return df.min(axis=1) < threshold


def get_kmer(pdb, threshold=0.45, scheme='closest-heavy'):
    """
    Computes the kmer sequence from a given PDB file based on contact mask.

    Parameters:
    pdb (str): The filename of the PDB file to be read.
    threshold (float, optional): The distance threshold to determine if a contact exists. Default is 0.45.
    scheme (str, optional): The scheme used for computing contacts. This could be any valid scheme recognized by the md.compute_contacts function from the MDTraj library. Default is 'closest-heavy'.

    Returns:
    str: A string representing the kmer sequence where contacts are indicated by the actual DNA base and non-contacts by '_'.
    """
    dna1 = get_dna1(pdb)
    mask = get_mask(pdb, threshold, scheme)
    return ''.join([x if y else '_' for x,y in zip(dna1, mask)])


def get_dna1(pdb):
    """
    Extracts the DNA sequence from a given PDB file name.

    Parameters:
    pdb (str): The filename of the PDB file.

    Returns:
    str: The extracted DNA sequence.
    """
    return pdb.split('---')[1].split('_')[2]
