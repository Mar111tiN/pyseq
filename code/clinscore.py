from yaml import CLoader as Loader, load
import pandas as pd
import numpy as np
from script_utils import show_output
from multiprocessing import Pool
from functools import partial

def load_weights(clinscore_file):
    '''
    load the relevant clinscore files into location and type dictionary
    '''
    
    with open(clinscore_file, "r") as stream:
        cosmic_weights = load(stream, Loader=Loader)
    return cosmic_weights


def cosmic_score_proc(df, clinscore_weights={}, verbose=2):
    '''
    computes the clinscore from a clinscore YAML file
    '''
    
    type_weights, loc_weights = clinscore_weights['type'], clinscore_weights['location']
    if verbose > 1:
        show_output(
            f"Computing cosmic score for {len(df.index)} mutations", multi=True)

    def compute_cosmic_score(row):
        """
        row-wise computation of cosmic90 scores
        """
        cos_score = 1 + type_weights.get(row["types"], 0) * loc_weights.get(row["location"], 0)
    
        return cos_score* int(row["count"])
    
    cosmic90_pattern = (r"(?P<count>[0-9]+)x\((?P<types>[^@)]+)(?:@(?P<location>[^0-9@)]+))?\)")
    # df = df['type'].str.replace("_(sclerosing_haemangioma)", "", regex=False).str.extractall(cosmic90_pattern)
    df["cosmic_score"]= df['type'].str.replace("_(sclerosing_haemangioma)", "", regex=False).str.extractall(cosmic90_pattern).apply(compute_cosmic_score, axis=1).reset_index().drop(columns="match").groupby("level_0").sum().fillna(0).astype(int)
    if verbose > 1:
        show_output("Finished", multi=True)
    return df


def get_cosmic_score(cosmic_df, cosmic_weights_file="", threads=2, verbose=1):
    '''
    computes the clinscore file for a "type" column using multicore processing
    '''

    # load the weights
    w= load_weights(cosmic_weights_file)
    # split the data 
    cosmic_split = np.array_split(cosmic_df, threads)
    if verbose:
        show_output(f"Computing cosmic score using {threads} threads.")
    pool = Pool(threads)
    dfs = pool.map(partial(cosmic_score_proc, clinscore_weights=w, verbose=verbose), cosmic_split)
    df = pd.concat(dfs)

    #### DEBUG
    # print(df.query("cosmic_score != cosmic_score"))
    df.loc[:, "cosmic_score"] = df["cosmic_score"].fillna(0).astype(int)

    # include the gene-wise multiplication factors if applicable
    gene_col = [col for col in df.columns if col.startswith("Gene")][0]

    if "genes" in w.keys():
        show_output(f'Inflating gene-wise scores for the following genes:\n{",".join(w["genes"].keys())}')
        for gene in w['genes']:
            df.loc[df[gene_col] == gene, "cosmic_score"] = df['cosmic_score'] * w['genes'][gene]
    
    if verbose:
        show_output(f"Cosmic score finished.")
    return df


def condense_muts_proc(df):
    '''
    takes a data frame with 'Chr', 'Start', 'End', 'Ref', 'Alt' and cosmic type and score and 
    removes duplicate mutations and aggregates different mutations at same positions
    other columns are just reduced to the first occurrence (might not be what you want), so..
    ... to save what you need, implement a sorting that keeps the keepers on top
    '''
    # first, remove duplicate exact mutations
    df = df.drop_duplicates(['Chr', 'Start', 'End', 'Ref', 'Alt'])   
    org_cols = df.columns
    # then, in one big swoop, keep the positions and the keepers of the other data columns
    # ..and combine with the aggregated data per mutations
    
    df = df.groupby(['Chr', 'Start', 'End']).first().reset_index().drop(['Ref', 'Alt', 'type', 'cosmic_score'], axis=1).merge(df.groupby(['Chr', 'Start', 'End']).agg({
    'Ref': lambda x: "/".join(x),
    'Alt': lambda x: "/".join(x),
    'type': lambda x: "+".join(x),
    'cosmic_score': 'sum' 
    }).reset_index()).loc[:, org_cols]  # restore original columns
    return df


def condense_muts_clinscore(df, threads, verbose=1):
    '''
    condenses the mutations to single hits per coord using condense_muts_proc
    '''

    if threads > 1:
        cosmic_split = np.array_split(df, threads)
        if verbose:
            show_output(f"Condensing cosmic scores using {threads} threads.")
        pool = Pool(threads)
        dfs = pool.map(condense_muts_proc, cosmic_split)    
        df = pd.concat(dfs)

    else:
        show_output(f"Condensing cosmic scores  on a single thread.")
        df = condense_muts_proc(df)
    show_output(f"Finished condensing cosmic scores.", color="success")
    return df