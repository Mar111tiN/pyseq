import pandas as pd
from script_utils import show_output
from pyseq_utils import full_collapse, remove_gene_dups, pos2bed
from clinscore import get_cosmic_score
from clinscore import condense_muts_clinscore


def roll(chrom_df, window_size):
    '''
    performs a rolling window computation for cosmic_score density and adds the results to the df
    '''
    
    # perform the rolling computation
    # Chr, Start, End are kept for computation as well as re-merging
    roll = chrom_df.rolling(window_size).agg({"Chr": "mean", "Start": min, "End":max, "cosmic_score": sum}).fillna(0).astype(int)
    # compute the density
    roll.loc[:, "cosmic_density"] = roll['cosmic_score'] / (roll['End'] - roll['Start'])
    # remerge (and drop preexising cosmic_density columns if any)
    chrom_df = chrom_df.drop(["cosmic_density"], axis=1, errors="ignore").merge(roll.drop(["End", "cosmic_score"], axis=1), on=['Chr', 'Start'], how="left")
    # fill NA overhangs via ffill
    chrom_df.loc[:, "cosmic_density"] = chrom_df["cosmic_density"].fillna(method="ffill")
    return chrom_df


def compute_cosmic_density(df, filter_setting={}, verbose=1):
    '''
    cycles through chromosomes and returns df with cosmic_density
    '''
    chrom_dfs = []
    if verbose:
        show_output("Computing mutation density")
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


def cosmic_panel_master(cosmic_df, cosmic_weights_file="", filter_setting={}, threads=10, verbose=1, condense_mut_positions=True):
    '''
    takes an annovar annotated mutation list and returns the collapsed mutation list based on filter list
    '''
    
    filter_info = "".join([f"\n\t[{col}:\t{filter_setting[col]}]" for col in ["cosmic_rolling_min", "rolling_window_size", "cosmic_min", "cosmic_density_min", "padding"]])
    show_output(f"Creating custom panel based on limits set in filter settings.{filter_info}")
    if cosmic_weights_file:
        # remove_duplicate_positions has to be set because we are only interested in the highest interest positions
        cosmic_scored = get_cosmic_score(cosmic_df, cosmic_weights_file=cosmic_weights_file, threads=threads, verbose=1)
        # reduce to unique mutations (by summing up the clinscore) (needed for panel design)
        # + group by start position and keep the first
        if condense_mut_positions:
            show_output("Condensing the mutations per position.")
            cosmic_scored = condense_muts_clinscore(cosmic_scored, threads=threads)
    else:
        if 'cosmic_score' in cosmic_df.columns:
            show_output(f"Using precomputed cosmic scores! For recomputation, provide a cosmic weights file")
            cosmic_scored = cosmic_df
        else:
            show_output("No clinscore in df and no weights file to compute clinscores. Sorry - stopping here!", color="warning")
            return

    # + sort by cosmic_score    # + sort by cosmic_score
    cosmic_scored = cosmic_scored.sort_values(['Chr', 'Start', 'cosmic_score'], ascending=[True, True, False])
    
    # perform rolling window computation
    if verbose:
        show_output("Perform rolling window computation")
    cosmic_denscored = compute_cosmic_density(cosmic_scored, filter_setting=filter_setting, verbose=verbose)

    # filter based on cosmic scores
    if verbose:
        show_output("Filtering out background mutations")
    panel_mut_df = filter_cosmic(cosmic_denscored, filter_setting=filter_setting, verbose=verbose)

    # collapse the df
    if verbose:
        show_output("Collapsing the mutations to adjacency groups")
    panel_mut_df, panel_region_df = full_collapse(panel_mut_df, padding=filter_setting['padding'], verbose=verbose)
    # meaningfull output
    mutN = panel_region_df['mutN'].sum()
    kb_size = int(panel_region_df['stretch'].sum() / 1000)
    show_output(f"Finished! Library size = {kb_size}kb - {mutN} mutations included", color="success")
    return panel_mut_df, panel_region_df, cosmic_denscored


def read_gene_list(panel_excel, gene_col="GeneInPanel"):
    '''
    reads a list of Genes ("GeneInPanel") with optional notes column and (n = XX) counts
    '''
    
    # read the excel file
    genes = pd.read_excel(panel_excel).dropna(subset=[gene_col]).reset_index(drop=True)
    genes = genes.loc[~genes[gene_col].str.startswith("#"), :]
    
    genes.loc[:, ['Gene', 'countOtherPanels']] = genes[gene_col].str.extract(r"(?P<Gene>[A-Z0-9]+)(?: n ? = ?)?(?P<countOtherPanels>[0-9]+)?")

    genes['countOtherPanels'] = genes['countOtherPanels'].fillna(1).astype(int)
    genes['notes'] = genes['notes'].str.strip(r"[ ()]").fillna("")
    genes = genes.sort_values(['Gene', 'notes']).groupby('Gene').agg({'notes':lambda x: "|".join(x), 'countOtherPanels':'sum'}).reset_index()
    genes['notes'] = genes['notes'].str.strip("|").str.replace(r"^[|]+$", "", regex=True)
    # remove possible Gene;Gene occurrences
    genes = remove_gene_dups(genes, gene_col='Gene')
    return genes.loc[:, ['Gene', 'countOtherPanels', 'notes']]


def analyze_genes(panel_mut_df, cosmic_scored, panel_excel="", save_excel=""):
    '''
    accumulate infos
    '''
    # get the top genes of all of cosmic
    top_genes = cosmic_scored.groupby("Gene").agg({'cosmic_score':"sum", 'type':'count'}).rename({'type':'count'}, axis=1).reset_index().sort_values('cosmic_score', ascending=False)
    
    # get the genes in the designed panel sorted by cosmic score
    gene_df = panel_mut_df.groupby("Gene").agg({'cosmic_score':'sum',  'type':'count'}).rename({'type':'count'}, axis=1).reset_index().sort_values('count', ascending=False)


    # merge top cosmic genes with genes from the designed panel
    merge = top_genes.merge(gene_df, on="Gene", how="outer", suffixes=('DB', 'Panel'))
    # edit columns
    for c in [c for c in merge.columns if "count" in c or "score" in c]:
        merge.loc[:, c] = merge[c].fillna(0).astype(int)
    
    # load the gene list from the published genes (panel_list)
    panel_gene_df = read_gene_list(panel_excel)

    # merge genes from other panels
    merge2 = merge.merge(panel_gene_df, on="Gene", how="outer")
    for c in [c for c in merge2.columns if "count" in c or "score" in c]:
        merge2.loc[:, c] = merge2[c].fillna(0).astype(int)
    merge2.loc[:, 'notes'] = merge2['notes'].fillna("").astype(str)
    
    # get 
    otherPanel_not_included = merge2.query('countOtherPanels > 0 and cosmic_scorePanel == 0').drop(['countDB', 'cosmic_scorePanel', 'countPanel'], axis=1)
    cosmic_not_included = merge2.query('cosmic_scoreDB > 50000 and cosmic_scorePanel == 0').drop(['cosmic_scorePanel', 'countPanel'], axis=1)
    in_panel = merge2.query('countPanel > 0')
    
    if save_excel:
        show_output(f"Saving to excel file {save_excel}.")
        # write to excel
        with pd.ExcelWriter(save_excel, mode="w") as writer:
            panel_mut_df.to_excel(writer, sheet_name="AllMutationsInPanel", index=False)
            pos2bed(panel_mut_df).to_excel(writer, sheet_name="AllMutationsBed", index=False, header=None)
            in_panel.to_excel(writer, sheet_name="GenesInPanel", index=False)
            cosmic_not_included.to_excel(writer, sheet_name="missing_TopCosmic", index=False)
            otherPanel_not_included.to_excel(writer, sheet_name="missing_otherPanels", index=False)
    
    return in_panel, cosmic_not_included, otherPanel_not_included