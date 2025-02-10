import numpy as np
from biopandas.pdb import PandasPdb
import mdtraj as md
import itertools
from scipy.stats import binomtest
import pandas as pd
from pandas_parallel_apply import DataFrameParallel
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import warnings
from glob import glob
import pymol2

def cif2pdb(infile):
    #Converting AF cifs to pdbs
    with pymol2.PyMOL() as pymol:
        pymol.cmd.feedback("disable", "all", "everything")
        pymol.cmd.load(infile, 'model_0',)
        pymol.cmd.save(infile.replace('.cif', '.pdb'), selection='model_0')


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
    pdb_file.read_pdb(pdb)

    # Creating dict with indexes of chains
    tdict = dict(
        pdb_file.amino3to1(fillna="").reset_index(drop=True).groupby(by='chain_id').apply(lambda x: x.index.tolist(), include_groups=False))

    # getting pairs
    pairs_A_C = list(itertools.product(tdict['A'], tdict['C']))
    pairs_B_C = list(itertools.product(tdict['B'], tdict['C']))
    pairs_A_D = list(itertools.product(tdict['A'], tdict['D']))
    pairs_B_D = list(itertools.product(tdict['B'], tdict['D']))

    file = md.load(pdb)

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


def get_kmer(pdb, dna, threshold=0.45, scheme='closest-heavy'):
    """
    Computes the kmer sequence from a given PDB file based on contact mask.

    Parameters:
    pdb (str): The filename of the PDB file to be read.
    threshold (float, optional): The distance threshold to determine if a contact exists. Default is 0.45.
    scheme (str, optional): The scheme used for computing contacts. This could be any valid scheme recognized by the md.compute_contacts function from the MDTraj library. Default is 'closest-heavy'.

    Returns:
    str: A string representing the kmer sequence where contacts are indicated by the actual DNA base and non-contacts by '_'.
    """
    mask = get_mask(pdb, threshold, scheme)
    return ''.join([x if y else '_' for x,y in zip(dna, mask)])


def get_kmer_rosetta(df):
    
    dna = df['DNA']
    pdb = f'{path}/commit_6a1b3d_pdbs/' + df['RoseTTAFold2NA_pdb_file']
    
    kmer = get_kmer(pdb, dna, threshold=threshold, scheme=scheme)
    
    return kmer


def get_kmer_alphafold(df):
    
    dna = df['DNA']
    pdb = df['AlphaFold3_pdb_file']
    
    pdbs = [f'{path}/alphafold3/{pdb}/{pdb}_model_{x}.pdb' for x in '01234']
    kmers = [get_kmer(x, dna, threshold=threshold, scheme=scheme) for x in pdbs]
    return kmers


def score_rosetta(df):
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
    kmers = df['RoseTTAFold2NA_extracted_kmers']

    # Check the mode and perform the comparison based on the suffix of the mode. Reverse case should end with either _r or _1. The length of flanks is 5. Comparison only takes place in the region of DNA that potentially binds to TFs, without considering flanks or the distance between binding sites
    if df['Side_of_preferred_DNA'] == "right":
        return kmers[5:ln+5].count("_") > kmers[-ln-5:-5].count("_")
    else:
        return kmers[5:ln+5].count("_") < kmers[-ln-5:-5].count("_")
    

def score_alphafold(df):
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
    five_kmers = df['AlphaFold3_extracted_kmers']

    # Check the mode and perform the comparison based on the suffix of the mode. Reverse case should end with either _r or _1. The length of flanks is 5. Comparison only takes place in the region of DNA that potentially binds to TFs, without considering flanks or the distance between binding sites
    if df['Side_of_preferred_DNA'] == "right":
        return [kmers[5:ln+5].count("_") > kmers[-ln-5:-5].count("_") for kmers in five_kmers]
    else:
        return [kmers[5:ln+5].count("_") < kmers[-ln-5:-5].count("_") for kmers in five_kmers]
    

def pvalues(df_input, side, column, rename):
    
    df_p = pd.crosstab(
        df_input[df_input['Side_of_preferred_DNA'] == side]['mode'], 
        df_input[df_input['Side_of_preferred_DNA'] == side][column],
        dropna=False)
    
    if True not in df_p.columns:
        df_p[True] = 0
    if False not in df_p.columns:
        df_p[False] = 0
    
    df_p = df_p.rename(columns={True : 'True', False : 'False'})[['True','False']]
    
    df_p = df_p.apply(lambda x : binomtest(x['True'], x['True'] + x['False']).pvalue, axis=1)
    
    df_p = df_p.reset_index().rename(columns={0 : rename})
    
    df_p['Side of the preferred DNA'] = side
    
    return df_p


def building_plot(df, df_f, dataset):

        # Sample sizes and significance levels (as provided earlier)
    sample_sizes = {"default": 119, "spac_m2": 96}

    # Define the figure and subplots
    fig, axes = plt.subplots(2, 1, figsize=(14, 12), sharex=True)

    # Custom names for scores and reverse descriptions
    score_labels = {'RoseTTAFold2NA_score': 'RoseTTAFold2NA', 'AlphaFold3_final_score': 'AlphaFold 3'}

    titles = {
        
        'left' : "Preferred DNA on the left side of 5'-3' input",
        'right' : "Preferred DNA on the right side of 5'-3' input"
    }

    # Plot for reverse False on top and True on the bottom
    for i, (reverse_value, title) in enumerate([('left', titles['left']), ('right', titles['right'])]):

        reverse_df = df[df["Side_of_preferred_DNA"] == reverse_value]
        
        # Separate data for RoseTTAFold2NA and AlphaFold 3
        rosetta_data = reverse_df['RoseTTAFold2NA_score']
        af_score_data = reverse_df['AlphaFold3_final_score']
        
        # Unique modes and corresponding positions
        modes = reverse_df['mode']
        x = np.arange(len(modes))
        bar_width = 0.4
        
        # Plot RoseTTAFold2NA and AlphaFold 3 side by side with transparency
        bars_rosetta = axes[i].bar(x - bar_width / 2, rosetta_data, width=bar_width, label=score_labels['RoseTTAFold2NA_score'], 
                                alpha=0.6, color='darkorange')
        bars_af = axes[i].bar(x + bar_width / 2, af_score_data, width=bar_width, label=score_labels['AlphaFold3_final_score'], 
                            alpha=0.6, color='orangered')
        
        # Titles and labels
        axes[i].grid(True)  # Enable grid
        axes[i].grid(color='lightgray', linestyle='-', linewidth=0.5)
        axes[i].set_ylabel("Fraction of Right Prediction", fontsize=14)
        axes[i].set_title(title, fontsize=14)
        axes[i].legend(loc='upper right')
        axes[i].set_ylim(0, 1)
        
        # Add a thicker line at y=0.5 with a label positioned above the line
        axes[i].axhline(0.5, color='gray', linestyle='--', linewidth=2)
        axes[i].text(-0.5, 0.52, '0.5', color='gray', ha='right', va='bottom', fontsize=10)
        axes[i].spines[['right', 'top']].set_visible(False)
        axes[i].tick_params(direction='in')


        # Annotate each bar with sample size (n) and significance asterisks only
        for bar, mode in zip(bars_rosetta, modes):
            
            n_value = sample_sizes["spac_m2"] if "spac_m2" in mode else sample_sizes["default"]
            
            p_value = df_f[(df_f['Side of the preferred DNA'] == reverse_value) & (df_f['mode'] == mode)]['RoseTTAFold2NA'].values[0]

            significance = '**' if p_value < 0.01 else '*' if p_value < 0.05 else ''
            
            # Adjust asterisk position for spac_p2 RoseTTAFold2NA
            asterisk_y_offset = 0.1 if (mode == 'spac_p2') else 0.08
            
            # Display n and asterisks (bold) labels
            axes[i].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + asterisk_y_offset, significance, 
                        ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')
            axes[i].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f'n={n_value}', 
                        ha='center', va='bottom', fontsize=10, color='black')
        
        for bar, mode in zip(bars_af, modes):
            
            n_value = sample_sizes["spac_m2"] if "spac_m2" in mode else sample_sizes["default"]
            
            p_value = df_f[(df_f['Side of the preferred DNA'] == reverse_value) & (df_f['mode'] == mode)]['AlphaFold 3'].values[0]
            
            significance = '**' if p_value < 0.01 else '*' if p_value < 0.05 else ''
            
            n_y_offset = 0.08 if (mode == 'orient_1_3' and reverse_value == 'right') else 0.02
            asterisk_y_offset_af = 0.1

            # Display n and asterisks (bold) labels
            axes[i].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + asterisk_y_offset_af, significance, 
                        ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')
            axes[i].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + n_y_offset, f'n={n_value}', 
                        ha='center', va='bottom', fontsize=10, color='black')

    # Set x-axis labels and layout adjustments
    axes[1].set_xlabel("Mode")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(modes)
    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(f'../results/supplementary_figure_{dataset}.svg', format='svg',dpi=300)


def building_plot_summary(df_summary, df_pvalues, dataset):
    
    # 1. Add overall geometry row (if not present)
    if 'overall geometry' not in df_summary['Side_of_preferred_DNA'].values:
        overall_rf = 6/7  # for RoseTTAFold2NA
        overall_af = 7/7  # for AlphaFold3_final
        
        df_overall = pd.DataFrame({
            'Side_of_preferred_DNA': ['overall geometry'],
            'RoseTTAFold2NA_score': [overall_rf],
            'AlphaFold3_final_score': [overall_af]
        })
        df_summary = pd.concat([df_summary, df_overall], ignore_index=True)
        
        # Compute overall p-values via binomtest (using p=0.5 as null hypothesis)
        pval_rf = binomtest(6, 7, p=0.5).pvalue
        pval_af = binomtest(7, 7, p=0.5).pvalue
        df_overall_p = pd.DataFrame({
            'Side_of_preferred_DNA': ['overall geometry'],
            'RoseTTAFold2NA_score': [pval_rf],
            'AlphaFold3_final_score': [pval_af]
        })
        df_pvalues = pd.concat([df_pvalues, df_overall_p], ignore_index=True)
    
    # 2. Reorder rows so that the order is:
    #    Overall geometry, Preferred DNA on the left, Preferred DNA on the right
    desired_order = ['overall geometry', 'left', 'right']
    df_summary = df_summary.set_index('Side_of_preferred_DNA').loc[desired_order].reset_index()
    df_pvalues = df_pvalues.set_index('Side_of_preferred_DNA').loc[desired_order].reset_index()
    
    # 3. Set up the bar plot parameters.
    # Two groups: "RF" and "AF"
    methods = ['RoseTTAFold2NA', 'AlphaFold 3']
    x = np.arange(len(methods))  # positions: 0 and 1
    
    n_subgroups = 3   # overall geometry, left, right
    bar_width = 0.2
    # Offsets so the three bars are centered at each x position.
    offsets = np.linspace(-bar_width, bar_width, n_subgroups)
    
    # Extract scores (in the desired order)
    scores_RF = df_summary['RoseTTAFold2NA_score'].values  # shape (3,)
    scores_AF = df_summary['AlphaFold3_final_score'].values   # shape (3,)
    
    # Extract p-values (in the same order)
    pvals_RF = df_pvalues['RoseTTAFold2NA_score'].values
    pvals_AF = df_pvalues['AlphaFold3_final_score'].values
    
    # Colors for each subgroup.
    # Original order was: [light, medium, dark].
    # Now we swap the light and dark so that overall geometry gets dark.
    colors_rf = ['#cc0000', '#ff6666', '#ffcccc']   # dark, medium, light red
    colors_af = ['#0033cc', '#6699ff', '#cce6ff']     # dark, medium, light blue
    
    # Define custom legend labels reflecting the new color mapping.
    legend_labels = [
        "Dark: Overall geometry",
        "Medium: Preferred DNA on the left",
        "Light: Preferred DNA on the right"
    ]
    
    # Helper function for significance stars.
    def significance_marker(p):
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return ''
    
    # 4. Create the plot.
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Loop over each subgroup (overall geometry, left, right)
    for i in range(n_subgroups):
        # RF bar: at group 0
        pos_RF = x[0] + offsets[i]
        ax.bar(pos_RF, scores_RF[i], bar_width, color=colors_rf[i])
        
        # AF bar: at group 1
        pos_AF = x[1] + offsets[i]
        ax.bar(pos_AF, scores_AF[i], bar_width, color=colors_af[i])
        
        # Annotate significance stars for both RF and AF bars.
        marker_RF = significance_marker(pvals_RF[i])
        if marker_RF:
            ax.text(pos_RF, scores_RF[i] * 1.01, marker_RF, ha='center', va='bottom', fontsize=14)
        marker_AF = significance_marker(pvals_AF[i])
        if marker_AF:
            ax.text(pos_AF, scores_AF[i] * 1.01, marker_AF, ha='center', va='bottom', fontsize=14)
    
    # Set x-axis ticks and labels.
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=12)
    
    # Apply a "Nature–like" style: remove grid and top/right spines.
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Set y-axis label and plot title.
    ax.set_ylabel('Fraction of the correct prediction', fontsize=12)
    ax.set_title('Prediction Performance: Overall Geometry & Spacing/Orientation Prediction', fontsize=14)
    
    # 5. Add the grey dashed horizontal line at y=0.5 and label it "0.5".
    ax.axhline(y=0.5, color='grey', linestyle='--', linewidth=1)
    x_min, x_max = ax.get_xlim()
    ax.text(x_min + 0.01*(x_max - x_min), 0.5, '0.5', color='grey',
            va='bottom', ha='left', fontsize=10)
    
    # 6. Create a custom legend with neutral patches.
    # Use neutral colors that mimic the new mapping: dark, medium, light.
    neutral_colors = ['#777777', '#999999', '#bbbbbb']
    custom_handles = [
        Patch(facecolor=neutral_colors[i], label=legend_labels[i])
        for i in range(n_subgroups)
    ]
    ax.legend(custom_handles, [h.get_label() for h in custom_handles],
              frameon=False, fontsize=12)
    
    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(f'../results/summary_figure_{dataset}.svg', format='svg',dpi=300)
    

def run(ncores, threshold_local, scheme_local, dataset):

    #Assigning dataset and global variables because DataFrameParallel has troubles working with lambda functions
    global path
    if dataset == 'sample':
        path = '../data/sample'
    elif dataset == 'full':
        path = '../data'
    global threshold 
    threshold = threshold_local
    global scheme 
    scheme = scheme_local

    #Converting AF3 cif files to pdb
    print('Converting CIF files of AlphaFold 3 to PDB')
    [cif2pdb(x) for x in sorted(glob(f'{path}/alphafold3/*/*.cif'))]

    #Reading DataFrame
    df_input = pd.read_csv(f'{path}/df_input.csv')

    #Extracting kmers from PDB files
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        print('Extracting kmers from AlphaFold 3 PDB files')
        df_input['AlphaFold3_extracted_kmers'] = DataFrameParallel(df_input, n_cores=ncores, pbar=False).apply(get_kmer_alphafold, axis=1)
        print('Extracting kmers from RoseTTAFold2NA PDB files')
        df_input['RoseTTAFold2NA_extracted_kmers'] = DataFrameParallel(df_input, n_cores=ncores, pbar=False).apply(get_kmer_rosetta, axis=1)

    #Calculating scores
    print('Calculating scores')
    df_input['RoseTTAFold2NA_score'] = df_input.apply(score_rosetta, axis=1)
    df_input['AlphaFold3_scores'] = df_input.apply(score_alphafold, axis=1)
    df_input['AlphaFold3_final_score'] = df_input['AlphaFold3_scores'].apply(lambda x : np.mean(x) > 0.5)
    df_input.to_csv(f'../results/df_input_filled_{dataset}.csv', index=False)

    #Calculating pvalues
    print('Calculating p values from binomial test')
    df_p1 = pvalues(df_input, 'left','RoseTTAFold2NA_score','RoseTTAFold2NA')
    df_p2 = pvalues(df_input, 'right','RoseTTAFold2NA_score','RoseTTAFold2NA')
    df_p3 = pvalues(df_input, 'left','AlphaFold3_final_score','AlphaFold 3')
    df_p4 = pvalues(df_input, 'right','AlphaFold3_final_score','AlphaFold 3')
    df_c1 = pd.concat([df_p1,df_p2],axis=0).reset_index(drop=True)
    df_c2 = pd.concat([df_p3,df_p4],axis=0).reset_index(drop=True)
    df_f = df_c1.merge(right=df_c2, on=['mode','Side of the preferred DNA'])
    df_f = df_f.iloc[:,[0,1,3,2]]
    df_f.to_csv(f'../results/pvalues_binom_{dataset}.csv', index=False)

    #Building plot
    print('Building bar chart with fractions of correct predictions')
    df_fractions = df_input.groupby(['Side_of_preferred_DNA','mode']).agg({
    'RoseTTAFold2NA_score': 'mean',  
    'AlphaFold3_final_score': 'mean'
    }).reset_index()
    building_plot(df=df_fractions, df_f=df_f, dataset=dataset)

    #Calculating summary
    #Calculating summary fractions
    summary = df_input.groupby('Side_of_preferred_DNA').agg({
    'RoseTTAFold2NA_score': 'mean',  
    'AlphaFold3_final_score': 'mean'
    }).reset_index()
    summary.to_csv(f'../results/summary_fractions_{dataset}.csv', index=False)
    #Calculating summary pvalues
    df_ps = df_input.groupby('Side_of_preferred_DNA').agg({
    'RoseTTAFold2NA_score': lambda x : binomtest(x.sum(),len(x)).pvalue,  
    'AlphaFold3_final_score': lambda x : binomtest(x.sum(),len(x)).pvalue
    }).reset_index()
    df_ps.to_csv(f'../results/summary_pvalues_{dataset}.csv', index=False)
    building_plot_summary(df_summary=summary, df_pvalues=df_ps, dataset=dataset)
    print('Summary was created')
