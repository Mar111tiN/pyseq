import pandas as pd


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


def filter_exonic(anno_df, filter_settings):
    '''
    filters the annotation list for relevant mutations
    '''
    
    ini_len = len(anno_df.index)
    # filter for exonic mutations
    exonic = anno_df['Func'].isin(filter_settings['exonic_list'])
    # filter for functional mutations
    SNV = anno_df['ExonicFunc'].isin(filter_settings['mut_list'])
    SNP = anno_df['gnomAD'].astype(float) > filter_settings['gnomad_max']
    cosmic_filtered = anno_df[exonic & SNV & ~SNP]
    filter_len = len(cosmic_filtered.index)
    print(f"Filtered out {ini_len - filter_len} mutations [{ini_len} --> {filter_len}]")
    return cosmic_filtered