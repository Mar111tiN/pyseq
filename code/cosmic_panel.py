import pandas as pd
from script_utils import show_output


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


def filter_cosmic(df, filter_settings={}, verbose=1):
    '''
    filters the mutation list based on cosmic score and mutation density
    '''
    ini_len = len(df.index)
    
    cosmic_min = filter_settings['cosmic_min']
    density_min = filter_settings['cosmic_density_min']
    df = df.query('(cosmic_score > @cosmic_min) or (cosmic_density > @density_min)')
    filter_len = len(df.index)
    if verbose:
        show_output(f"Filtered out {ini_len - filter_len} mutations [{ini_len} --> {filter_len}]")
    return df


def collapse(chrom_df, pad=100):
    '''
    detect overlapping regions and collapse stretches
    expects single chromosome dfs
    '''
    pd.options.mode.chained_assignment = None
    # copy just to make sure
    cr = chrom_df.copy()
    #extend the coords
    cr['Start'] = cr['Start'] - pad
    cr['End'] = cr['End'] + pad
    # get the overlaps
    cr['ov1'] = (cr['End'] > cr.shift(-1)['Start']).astype(int)
    cr['ov2'] = (cr['Start'] < cr.shift(1)['End']).astype(int)
    # assign overlap groups
    cr['ovgroup'] = ((cr['ov1'] * (cr['ov2'] == 0).astype(int)).cumsum()) * (cr['ov1'] | (cr['ov2']))   
    # remove the isolated mutations from the grouping because they are all group 0 and re-index
    cr_nogap = cr.query('ovgroup == 0')
    cr = cr.query('ovgroup > 0').reset_index(drop=True)
    # assign negative indices to isolated mutations
    cr_nogap.loc[:,'ovgroup'] = -cr_nogap.index
    # re-combine for proper grouping
    cr = pd.concat([cr_nogap, cr]).reset_index(drop=True).sort_values(["Chr", "Start"], ascending=True)
    # condense the groups and keep important metrices
    cg = cr.groupby("ovgroup").agg(dict(
        Chr="first",
        Start="min",
        End="max",
        Gene=["first","last"],
        cytoband="min",
        gnomAD="max",
        cosmic_score="sum",
        cosmic_density="mean",
        Func="count"
    ))
    
    # reassign the multiindex to simple column index
    cg.columns = [f"{col[0]}-{col[1]}" if col[0] == "Gene" else col[0] for col in cg.columns]
    # consolidate the Gene Info
    cg.loc[cg['Gene-first'] == cg['Gene-last'], 'Gene-last'] = ""
    cg = cg.rename({"Gene-first": "Gene", "Gene-last":"Gene2", "Func":"mutN"}, axis=1).reset_index()
    # add the length of the overlap
    cg.loc[:, 'stretch'] = cg['End'] - cg['Start']
    #reorder columns
    # cg.columns = list(cg.columns[1:]) + [cg.columns[0]]
    cg = cg.sort_values(['Chr', 'Start'], ascending=True).reset_index(drop=True)
    pd.options.mode.chained_assignment = "warn"
    return cg, cr

def full_collapse(df, padding=100, verbose=1):
    '''
    cycles through chromosomes and returns df with collapsed mutation regions
    '''
    if verbose:
        show_output("Collapsing adjacent mutations and including bait padding")
    chrom_dfs = []
    chrom_group_dfs = []
    for chrom in df['Chr'].unique():
        if verbose > 1:
            show_output(f"Collapsing chromosome {chrom}", time=False)
        chrom_group_df, chrom_df = collapse(df.query("Chr == @chrom"), pad=padding)
        chrom_group_dfs.append(chrom_group_df)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs)
    group_df = pd.concat(chrom_group_dfs).loc[:,['Chr', 'Start', 'End', 'Gene', 'Gene2', 'cytoband', 'gnomAD',
       'cosmic_score', 'cosmic_density', 'ovgroup', 'mutN', 'stretch']].sort_values(['Chr', 'Start'])
    return df, group_df