import pandas as pd
import os
from script_utils import show_output


def remove_gene_dups(df, gene_col="Gene"):
    '''
    cosmic dbs sometimes contain Gene entries of this kind: RB1;RB1
    this tool removes these entries
    '''
    for i in range(3):
        df.loc[df[gene_col].str.contains(";"), gene_col] = df[gene_col].str.replace(r"([A-Z0-9]+);\1", r"\1", regex=True)
    return df


def remove_cosmic_dups(df):
    '''
    for some reason, sometimes genes appear more often in the DB with different types
    probably after applying remove_gene_dups
    '''
    no_dups = df.loc[~df.duplicated(['Chr', 'Start', 'End', 'Ref', 'Alt'], keep=False), :]
    dups = df.loc[df.duplicated(['Chr', 'Start', 'End', 'Ref', 'Alt'], keep=False), :]
    # remove by grouping and concating the types
    dups = dups.groupby([c for c in dups.columns if c not in ["type", "Mut_ID"]]).agg({"Mut_ID": "first", "type":lambda x: "+".join(x)}).reset_index()
    df = pd.concat([no_dups, dups]).sort_values(['Chr', 'Start', 'End', 'Ref', 'Alt'])
    return df



# load the annovar output
def load_anno(file):
    '''
    load the annovar file and edits columns
    '''
    cosmic_anno = pd.read_csv(file, sep="\t", low_memory=False)
    # quick edit of the columns
    cosmic_anno.columns = [cosmic_anno[col].iloc[0] if col.startswith("Other") else col.replace(".refGene","").replace("_exome_ALL", "") for col in cosmic_anno.columns]
    cosmic_anno = cosmic_anno[1:]
    # adjust file types
    for col in ["Chr", "Start", "End"]:
        cosmic_anno.loc[:,col] = cosmic_anno[col].fillna(0).astype(float).astype(int)
        cosmic_anno.loc[cosmic_anno['gnomAD'] == ".",'gnomAD'] = 0
    cosmic_anno.loc[:, 'gnomAD'] = cosmic_anno['gnomAD'].astype(float)
    # remove unneeded cols
    cosmic_anno = cosmic_anno.loc[:,['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func', 'Gene',
       'ExonicFunc', 'AAChange', 'cytoband', 'gnomAD', 'Mut_ID', 'type']]
    return cosmic_anno


def filter_exonic(anno_df, filter_setting):
    '''
    filters the annotation list for relevant mutations
    '''
    
    ini_len = len(anno_df.index)
    # filter for exonic mutations
    exonic = anno_df['Func'].isin(filter_setting['exonic_list'])
    # filter for functional mutations
    SNV = anno_df['ExonicFunc'].isin(filter_setting['mut_list'])
    SNP = anno_df['gnomAD'].astype(float) > filter_setting['gnomad_max']
    cosmic_filtered = anno_df[exonic & SNV & ~SNP]
    filter_len = len(cosmic_filtered.index)
    print(f"Filtered out {ini_len - filter_len} mutations [{ini_len} --> {filter_len}]")
    return cosmic_filtered



def collapse(chrom_df, pad=100):
    '''
    detect overlapping regions and collapse stretches
    expects single chromosome dfs
    '''
    pd.options.mode.chained_assignment = None
    # copy just to make sure
    cr = chrom_df.copy()
    #extend the coords
    cr.loc[:,'Start'] = cr['Start'] - pad
    cr.loc[:,'End'] = cr['End'] + pad
    # get the overlaps
    cr.loc[:,'ov1'] = (cr['End'] > cr.shift(-1)['Start']).astype(int)
    cr.loc[:,'ov2'] = (cr['Start'] < cr.shift(1)['End']).astype(int)
    # assign overlap groups
    cr.loc[:,'ovgroup'] = ((cr['ov1'] * (cr['ov2'] == 0).astype(int)).cumsum()) * (cr['ov1'] | (cr['ov2']))   
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
            show_output(f"Collapsing chromosome {chrom}")
        chrom_group_df, chrom_df = collapse(df.query("Chr == @chrom"), pad=padding)
        chrom_group_dfs.append(chrom_group_df)
        chrom_dfs.append(chrom_df)
    df = pd.concat(chrom_dfs).drop(['ov1', 'ov2'], axis=1)
    group_df = pd.concat(chrom_group_dfs).loc[:,['Chr', 'Start', 'End', 'Gene', 'Gene2', 'cytoband', 'gnomAD',
       'cosmic_score', 'cosmic_density', 'ovgroup', 'mutN', 'stretch']].sort_values(['Chr', 'Start'])

    # now, remove the padding again
    df.loc[:,'Start'] = df['Start'] + padding
    df.loc[:,'End'] = df['End'] - padding
    return df, group_df


def collapse_bed(chrom_df):
    '''
    detect overlapping regions and collapse stretches
    expects single chromosome dfs
    '''
    pd.options.mode.chained_assignment = None
    # copy just to make sure
    cr = chrom_df.copy()
    #extend the coords
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
    ))
    
    # consolidate the Gene Info
    # add the length of the overlap
    cg.loc[:, 'stretch'] = cg['End'] - cg['Start']
    #reorder columns
    # cg.columns = list(cg.columns[1:]) + [cg.columns[0]]
    cg = cg.sort_values(['Chr', 'Start'], ascending=True).reset_index()
    pd.options.mode.chained_assignment = "warn"
    return cg


def is_bed(bed_df):
    '''
    checks for chromosome label and Start-End id
    '''
    
    test = bed_df.iloc[0,:]
    try:
        has_chr = test[0].startswith("chr") or int(test[0]) < 25
    except:
        has_chr = False      
    try:
        start_end = test[1].is_integer() and test[2].is_integer()
    except:
        start_end = False                                  
    return has_chr and start_end


def get_bedsize(bed_file, verbose=0):
    '''
    reads bedfile and returns the library size
    cycles through chromosomes and returns df with collapsed mutation regions
    '''
    # read bed_file and skip improper rows
    i = 0
    while True:
        df = pd.read_csv(bed_file, sep="\t", skiprows=i, names=['Chr', 'Start', 'End', 'ID', 'desc', 'info'])
        if is_bed(df):
            show_output(f"Skipping first {i} rows of bedfile", color="warning")
            break
        i += 1
        
        
    if verbose:
        show_output("Collapsing adjacent mutations")
    chrom_group_dfs = []
    for chrom in df['Chr'].unique():
        if verbose > 1:
            show_output(f"Collapsing chromosome {chrom}")
        chrom_group_df = collapse_bed(df.query("Chr == @chrom"))
        chrom_group_dfs.append(chrom_group_df)
    group_df = pd.concat(chrom_group_dfs).loc[:,['Chr', 'Start', 'End', 'ovgroup','stretch']].sort_values(['Chr', 'Start'])
    
    bedsize = group_df['stretch'].sum()
    show_output(f"bedfile {os.path.basename(bed_file)} has a design size of {bedsize / 1000}kb", color="success")
    return bedsize


def pos2bed(df, chr_start_end=[], as_string=False):
    '''
    takes a dataframe with three columns Chr, Start and End and converts and returns a 
    pd.Series() containing the respective bed coords as "Chr:Start-End"
    If the input file has different cols for Chr, Start, End provide them as list of col names
    '''

    if len(chr_start_end) == 3:
        coords = df[chr_start_end[0]] + ":" + (df[chr_start_end[1]]).astype(str) + "-" + (df[chr_start_end[2]]).astype(str)
    else:
        coords = df['Chr'] + ":" + (df['Start']).astype(str) + "-" + (df['End']).astype(str)
    # convert to string without indices for direct copy to clipboard
    if as_string:
        coords = coords.to_string(index=False)
    return coords


def bed2pos(col):
    '''
    converts a df column into a df with the three columns Chr Start End
    '''
    s = col.str.extract(r"(?P<Chr>chr[0-9]+):(?P<Start>[0-9]+)-(?P<End>[0-9]+)")

    return s