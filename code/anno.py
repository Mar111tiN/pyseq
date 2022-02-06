import os
import pandas as pd
import shutil
from yaml import CLoader as Loader, load

from script_utils import run_cmd, show_output, show_cmd


def get_anno_params(config_file, threads=4):
    '''
    helper function to create full annovar parameters from configs and threads
    the actual file will be included by replacing the <file> placeholder
    '''

    # load annovar configs
    with open(config_file, "r") as stream:
        config = load(stream, Loader=Loader)
    config
    
    # get the full path to humandb
    humandb = os.path.join(os.environ['STATIC'], config['humandb'])
    
    # get the available anno files
    file_list = [f for f in list(os.walk(humandb))[0][-1] if f.endswith(".txt")]
    # reduce anno files to the files compatible with genome build version
    build = config['build']
    build_files = []
    for file in file_list:
        if file.startswith(build):
            build_files.append(file)

    # filter the anno protocol for the available files for that genome build        
    anno_refs = config['annofiles']
    anno_list = []
    missing_list = []
    for anno in anno_refs:
        for file in build_files:
            if anno in file:
                anno_list.append(anno)
                break
        # if anno has not been found in file during for-loop
        else:
            missing_list.append(anno)

    # create the protocol string
    protocol = ','.join(anno_list)
    if missing_list:
        show_output(f"{' '.join(missing_list)} not found for {build}! Doing without.. ", color="warning")
    # create the operation string 'g,r,f,f,f,f' assuming all but the first three dbs (ref, cytoBand, superDups) in config to be filter-based
    operation_list = []
    for anno in anno_list:
        if "Gene" in anno:
            operation_list.append('g')
        elif anno in ['cytoBand', 'genomicSuperDups']:
            operation_list.append('r')
        else:
            operation_list.append('f')
    operation = ','.join(operation_list)
    
    # now build the command
    c_anno = f"perl {config['annovar_path']}/table_annovar.pl"
    c_fix = '-nastring "." --remove'
    c_thread = f"--maxgenethread {threads} --thread {threads}"
    full_cmd = f'{c_anno} -buildver {build} {c_thread} --protocol {protocol} --operation {operation} {c_fix} --outfile <out_file> <in_file> {humandb}'
    return full_cmd


def run_annovar(file, annovar_config={}, threads=6, temp_folder="", cleanup=True):
    '''
    runs the annovar command from a tab_separated file or a df using an annovar_config yaml and returning a df with the annotations attached
    if temp_folder is not provided, a temp folder is created (and deleted) 
        - either at the basedir of the file
        - or at the execution dir if df was provided
    data frame or file is expected to have Chr, Start, End, Ref, Alt columns
    '''
    
    is_df = isinstance(file, pd.DataFrame)

    ## create temp folder and load file or df
    if not temp_folder:
        if is_df:
            temp_folder = os.path.join(os.getcwd(), "anno_temp")    
        else:
            temp_folder = os.path.join(os.path.dirname(file), "anno_temp")
    # create the temp folder
    if not os.path.exists(temp_folder):
        try:
            os.makedirs(temp_folder)
        except:
            show_output(f"Temp folder {temp_folder} could not be created. Please check permissions!", color="warning")
            exit()
            
    # store file as df to better handle headers and extra columns
    df = file if is_df else pd.read_csv(file, sep="\t")
    for col in ["Start", "End"]:
        df.loc[:, col] = df[col].astype(int)
    # keep other columns in other_df
    other_cols = [col for col in df.columns if not col in ["Chr", "Start", "End", "Ref", "Alt"]]
    
    
    # store headerless file for annovar computation
    file = os.path.join(temp_folder, "anno_file.csv")
    df.loc[:,["Chr", "Start", "End", "Ref", "Alt"]].to_csv(file, sep="\t", index=False, header=False)
    
    out_file = os.path.join(temp_folder, "annovar_out")
    
    # remove header because otherwise annovar will use header as first row
    # i_nohead = f"{output_base}.nohead.csv"
    # run_cmd(f"mawk 'NR > 1 {{print}}' < {i} > {i_nohead}")

    anno_cmd = get_anno_params(annovar_config, threads=threads).replace("<out_file>", out_file).replace("<in_file>", file)
    run_cmd(anno_cmd, multi=False)

    # reload the file 
    anno_df = pd.read_csv(os.path.join(temp_folder, "annovar_out.hg38_multianno.txt"), sep="\t")
    # format_cmd = f"cat {output_base}.*.txt | {annoinfo} > {o}"
    # run_cmd(format_cmd)
    # remerge other_df
    if other_cols:
        anno_df = anno_df.merge(df, on=["Chr", "Start", "End", "Ref", "Alt"])
    
    # cleanup
    if cleanup:
        shutil.rmtree(temp_folder)
    return anno_df