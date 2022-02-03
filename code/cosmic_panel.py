import pandas as pd
from script_utils import show_output
from pyseq_utils import full_collapse
from clinscore import get_cosmic_score

def roll(chrom_df, window_size):
    '''
    performs a rolling window computation for cosmic_score density and adds the results to the df
    '''
    
    # perform the rolling computation
    # Chr, Start, End are kept for computation as well as re-merging
    roll = chrom_df.rolling(window_size).agg({"Chr": "mean", "Start": min, "End":max, "cosmic_score": sum}).fillna(0).astype(int)
    # compute the density
    roll.loc[:, "cosmic_density"] = roll['cosmic_score'] / (roll['End'] - roll['Start'])
    # remerge
    chrom_df = chrom_df.merge(roll.drop(["End", "cosmic_score"], axis=1), on=['Chr', 'Start'], how="left")
    # fill NA overhangs via ffill
    chrom_df.loc[:, "cosmic_density"] = chrom_df["cosmic_density"].fillna(method="ffill")
    return chrom_df


def compute_cosmic_density(df, filter_setting={}, verbose=1):
    '''
    cycles through chromosomes and returns df with cosmic_density
    '''
    chrom_dfs = []
    if verbose:
        show_output("Computing mutation density", time=False)
    # remove background mutations
    cosmin = filter_setting['cosmic_rolling_min']
    df = df.query('cosmic_score >= @cosmin')
    # cycle through chromosomes
    for chrom in df['Chr'].unique():
        if verbose > 1:
            show_output(f"Rolling on chromosome {chrom}")
        chrom_df = roll(df.query("Chr == @chrom"), window_size=filter_setting['rolling_window_size'])
        chrom_dfs.append(chrom_df)

    df = pd.concat(chrom_dfs)
    df.loc[:, 'cosmic_density'] = df['cosmic_density'].round(1)
    df.loc[:, "cosmic_score"] = df['cosmic_score'].astype(int)
    return df


def filter_cosmic(df, filter_setting={}, verbose=1):
    '''
    filters the mutation list based on cosmic score and mutation density
    '''
    ini_len = len(df.index)
    
    cosmic_min = filter_setting['cosmic_min']
    density_min = filter_setting['cosmic_density_min']
    df = df.query('(cosmic_score > @cosmic_min) or (cosmic_density > @density_min)')
    filter_len = len(df.index)
    if verbose:
        show_output(f"Filtered out {ini_len - filter_len} mutations [{ini_len} --> {filter_len}]")
    return df


def cosmic_master(df, cosmic_weights_file="", filter_setting={}, verbose=1):
    '''
    takes an annovar annotated mutation list and returns the collapsed mutation list based on filter list
    '''
    
    filter_info = "".join([f"\n\t[{col}:\t{filter_setting[col]}]" for col in ["cosmic_rolling_min", "rolling_window_size", "cosmic_min", "cosmic_density_min", "padding"]])
    show_output(f"Creating custom panel based on limits set in filter settings.{filter_info}")
    if cosmic_weights_file:
        df_scored = get_cosmic_score(df, cosmic_weights_file=cosmic_weights_file, threads=10, verbose=1)
    else:
        if 'cosmic_score' in df.columns:
            show_output(f"Using precomputed cosmic scores! For recomputation, provide a cosmic weights file", time=False)
            df_scored = df
        else:
            show_output("No clinscore in df and no weights file to compute clinscores. Sorry - stopping here!", color="warning")
            return
    
    # perform rolling window computation
    if verbose:
        show_output("Perform rolling window computation", time=False)
    df = compute_cosmic_density(df_scored, filter_setting=filter_setting, verbose=verbose)

    # filter based on cosmic scores
    if verbose:
        show_output("Filtering out background mutations", time=False)
    df = filter_cosmic(df, filter_setting=filter_setting, verbose=verbose)

    # collapse the df
    if verbose:
        show_output("Collapsing the mutations to adjacency groups", time=False)
    df, group_df = full_collapse(df, padding=filter_setting['padding'], verbose=verbose)

    # meaningfull output
    mutN = group_df['mutN'].sum()
    kb_size = int(group_df['stretch'].sum() / 1000)
    show_output(f"Library size = {kb_size}kb - {mutN} mutations included")
    return df, group_df, df_scored