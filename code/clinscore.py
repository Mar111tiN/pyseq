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
    
    # reduce to unique mutations (with highest clinscore)
    # + sort by cosmic_score
    # + group by start position and keep the first
    df = df.sort_values(['Chr', 'Start', 'cosmic_score'], ascending=[True, True, False])
    df = df.groupby(['Chr', 'Start']).first().reset_index()
    if verbose > 1:
        show_output("Finished", multi=True)
    return df


def get_cosmic_score(cosmic_df, cosmic_weights_file="", threads=2, verbose=1):
    '''
    computes the clinscore file for a "type" column using multicore processing
    '''

    # load the weights
    clinscore_weights = load_weights(cosmic_weights_file)
    # split the data 
    cosmic_split = np.array_split(cosmic_df, threads)
    if verbose:
        show_output(f"Computing cosmic score using {threads} threads.")
    pool = Pool(threads)
    cosmic_score_dfs = pool.map(partial(cosmic_score_proc, clinscore_weights=clinscore_weights, verbose=verbose), cosmic_split)
    cosmic_score_df = pd.concat(cosmic_score_dfs)

    #### DEBUG
    # print(cosmic_score_df.query("cosmic_score != cosmic_score"))
    cosmic_score_df.loc[:, "cosmic_score"] = cosmic_score_df["cosmic_score"].astype(int)
    if verbose:
        show_output(f"Cosmic score finished.")
    return cosmic_score_df