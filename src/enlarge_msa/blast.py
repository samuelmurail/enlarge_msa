import os 
import pandas as pd
import numpy as np
import subprocess
import glob
import shutil
from tqdm import tqdm
from multiprocessing.dummy import Pool
from concurrent.futures import ThreadPoolExecutor
from threading import Semaphore

def process_csv(csv_path):
    """
    Prepare a DataFrame from an Alphafold input CSV file.

    This function reads an Alphafold input CSV file to extract and prepare a DataFrame containing 
    separated peptide and receptor sequences.

    Parameters
    ----------
    csv_path : str
        The path to the Alphafold input CSV file.

    Returns
    -------
    df : pandas.DataFrame or None
        DataFrame containing columns 'pdb', 'rec_seq', and 'pep_seq' if the CSV file is processed successfully.
        None if there is an issue with the CSV file.
    """
    AF_df = pd.read_csv(csv_path)
    
    pep_seq = []
    unique_rec_seq = []
    full_rec_seq=[]
    pdb_list = AF_df["id"].tolist()
    
    for seq in AF_df["sequence"].tolist():
        pep = seq.split(":")[-1]
        rec = seq.split(":")[0:-1]
        if len(rec)!= len(set(rec)):
            full_rec_seq.append(rec)
        else : 
            full_rec_seq.append(None)
        pep_seq.append(pep)
        unique_rec_seq.append(pd.unique(np.array(rec)).tolist())
    
    if len(pep_seq) == len(unique_rec_seq) == len(pdb_list):
        data = {
            "pdb": pdb_list,
            "full_rec_seq":full_rec_seq,
            "unique_rec_seq": unique_rec_seq,
            "pep_seq": pep_seq,
        }
        
        df = pd.DataFrame(data)
        return df
    else:
        print(f"Error: problem with the csv. \n")
        return None
        
def blast_fasta_generator(rec_seq,pep_seq,pdb):

    """
    Generate FASTA files for receptor and peptide sequences.

    This function creates FASTA files for the given receptor and peptide sequences using the provided PDB .
    The files are saved in the current directory with filenames formatted as 'rec_<pdb>.fa' and 'pep_<pdb>.fa'.

    Parameters
    ----------
    rec_seq : str
        The receptor sequence to be saved in the FASTA file.
    pep_seq : str
        The peptide sequence to be saved in the FASTA file.
    pdb : str
        The PDB id of the query.

    Returns
    -------
    None
    """
    folder_name = './fasta_files'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    pep_file_path = os.path.join(folder_name, f"pep_{pdb}.fa")
    print(f"Preparing peptide FASTA file of {pdb} \n")
    
    with open(pep_file_path, "w") as file:
        file.write(f">pep_{pdb}\n")
        file.write(pep_seq)
        
    for i, seq in enumerate(rec_seq):
        rec_file_path = os.path.join(folder_name, f"rec_{pdb}_{i}.fa")
        print(f"Preparing receptor FASTA file of {pdb} \n")
        
        with open(rec_file_path, "w") as file:
            file.write(f">rec_{pdb}_{i}\n")
            file.write(seq)
            
def blasp_launch(fasta, db, output, num_threads, evalue,semaphore):
    """
    Launch a BLASTP search.

    This function runs a BLASTP search using the specified parameters. It constructs the BLASTP command line
    and executes it as a subprocess.

    Parameters
    ----------
    fasta : str
        The path to the input FASTA file containing the sequence to be queried.
    db : str
        The BLAST database to search against.
    output : str
        The path to the output file where the BLASTP results will be saved.
    num_threads : int
        The number of threads to use for the BLASTP search.
    evalue : float
        The e-value threshold for reporting matches.
    semaphore : threading.Semaphore
        A semaphore to limit the number of concurrent BLASTP processes.

    Returns
    -------
    None
    """
    with semaphore:
        BLASTP_BIN = shutil.which('blastp')
        if BLASTP_BIN is None:
            print(f'blastp could not be found \n')
            return

        command_line = [
            BLASTP_BIN, '-query', fasta, '-db', db, '-out', output,
            '-num_threads', str(num_threads), '-evalue', str(evalue) ,'-outfmt', '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps qcovs'
        ]
        subprocess.run(command_line)
        
def blastp(df, db, evalue, num_threads, pep=True, rec=False, max_runs=5):
    """
    Generate and run BLASTP searches for sequences in a DataFrame.

    This function generates FASTA files for receptor and peptide sequences from a DataFrame and runs BLASTP searches 
    on them using the specified database and parameters.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns 'rec_seq', 'pep_seq', and 'pdb'.
    db : str
        The path to the BLAST database to be searched against.
    evalue : float
        The e-value threshold for the BLASTP search.
    num_threads : int
        The number of threads to use for the BLASTP search.
    pep : bool, optional
        If True, perform BLASTP search on peptide sequences (default is True).
    rec : bool, optional
        If True, perform BLASTP search on receptor sequences (default is False).
    max_runs : int, optional
        The number of concurrent threads to use for the BLASTP searches (default is 5).

    Returns
    -------
    None
    """
    semaphore = Semaphore(max_runs)
    tasks_pep = []
    tasks_rec = []

    fasta_folder = './fasta_files'
    if not os.path.exists(fasta_folder):
        print(f"Fasta files folder is missing \n")

    blast_output = './blast_output'
    if not os.path.exists(blast_output):
        os.makedirs(blast_output)

    for rec_seq, pep_seq, pdb in (zip(df["unique_rec_seq"], df["pep_seq"], df["pdb"])):
        blast_fasta_generator(rec_seq, pep_seq, pdb)  
        if pep:
            print(f"Running BLASTp for peptide of: {pdb} \n")
            pep_file_path = os.path.join(fasta_folder, f"pep_{pdb}.fa")
            pep_blast_output=os.path.join(blast_output, f'out_pep_{pdb}.csv')
            tasks_pep.append((pep_file_path, db, pep_blast_output, num_threads, evalue, semaphore))

        if rec:
            print(f"Running BLASTp for receptors of :{pdb} \n")
            for i, seq in enumerate(rec_seq):
                rec_file_path = os.path.join(fasta_folder, f"rec_{pdb}_{i}.fa")
                rec_blast_output=os.path.join(blast_output, f'out_rec_{pdb}_{i}.csv')
                tasks_rec.append((rec_file_path, db, rec_blast_output, num_threads, evalue, semaphore))
                
    with ThreadPoolExecutor(max_workers=max_runs) as executor:
        if pep:
            executor.map(lambda p: blasp_launch(*p), tasks_pep)
        if rec:
            executor.map(lambda p: blasp_launch(*p), tasks_rec)

def blastdbcmd(query, db):
    """
    Fetch an amino acid sequence from the database.

    This function runs a BLAST+ blastdbcmd search to retrieve the amino acid sequence 
    for the specified query ID from the given BLAST database.

    Parameters
    ----------
    query : str
        The query ID.
    db : str
        The path to the BLAST database to be searched against.

    Returns
    -------
    str
        The amino acid sequence of the searched query, or an empty string if not found.
    """
    
    Blastdbcmd_BIN = shutil.which('blastdbcmd')
    if Blastdbcmd_BIN is None:
        print(f'blastdbcmd could not be found \n')
        return ""
        
    command_line = [Blastdbcmd_BIN, '-db', db, '-entry', query]
    result = subprocess.run(command_line, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running blastdbcmd: {result.stderr} \n")
        return ""
    
    lines = result.stdout.splitlines()
    seq = "".join(lines[1:]) if len(lines) > 1 else ""

    return seq
    
def process_blast_results(csv_path, seq, db):

    """
    Treat BLAST results and augment with sequences.

    This function reads a CSV file of BLAST search results and retrieves the corresponding 
    sequences from the BLAST searched database, and adds these sequences to the DataFrame.

    Parameters
    ----------
    csv_path : str
        The path to the CSV file containing BLAST search results.
    seq : str
        The query sequence used in the BLAST search.
    db : str
        The path to the BLAST database to retrieve sequences from.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the original BLAST search results, with additional column for 
        sequences.
    """
    
    blastp_df=None
    column_names = ["query", "subject_id", "% identity", "alignment_length", "mismatches", 
                "gap_open", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score" ,"query_len" , "subj_len" , "gaps" , "qcovs"]
    seq_list = []
    blastp_df = pd.read_csv(csv_path, names=column_names)        
    for query in blastp_df["subject_id"].tolist():
        new_seq = blastdbcmd(query, db)
        if new_seq:
            seq_list.append(new_seq)
        else:
            seq_list.append("") 
    
    blastp_df["sequence"] = seq_list    
    return blastp_df

def best_seq(blastp_df):
    """
    Identify and return the best sequence based on maximum coverage and sequence length.

    This function finds first the sequence with the highest coverage percentage in the provided DataFrame.
    If there are multiple sequences with the same highest coverage, it returns the longest sequence.

    Parameters
    ----------
    blastp_df : pd.DataFrame
        The DataFrame containing BLAST search results with additional column for sequence information.

    Returns
    -------
    str
        The sequence with the highest coverage percentage and the longest length. Returns an empty 
        string if the DataFrame is empty.
    """
    
    if blastp_df.empty:
        return ""

    max_coverage = blastp_df['qcovs'].max()
    max_coverage_df = blastp_df[blastp_df['qcovs'] == max_coverage]
    sequence_with_max_coverage_and_length_df = blastp_df.loc[max_coverage_df['subj_len'].idxmax()]
    new_seq = sequence_with_max_coverage_and_length_df["sequence"]

    return new_seq

def reconstruct(full_rec_seq, unique_rec_seq, new_rec_seq):
    """
    Reconstruction of a sequence.

    This function helps in reconstructing an original sequence with new best retrieved ones 
    for the case of proteins with repetitive sequences. It replaces each segment in the 
    original sequence with new ones, in order to obtain the final full sequence.

    Parameters
    ----------
    full_rec_seq : list of str
        The original sequence containing repetitive segments to be replaced.
    unique_rec_seq : list of str
        A list of unique segments obtained from the original sequence.
    new_rec_seq : list of str
        A list of new sequences that will replace the corresponding segments in the original one.

    Returns
    -------
    list of str
        The reconstructed sequence with segments replaced as specified.
        
    Examples
    --------
    >>> full_rec_seq = ['A', 'B', 'A']
    >>> unique_rec_seq = ['A', 'B']
    >>> new_rec_seq = ['X', 'Y']
    >>> reconstruct(full_rec_seq, unique_rec_seq, new_rec_seq)
    ['X', 'Y', 'X']
    """

    reconstructed_seq = []
    for seq in full_rec_seq:
        for unique_seq, new_seq in zip(unique_rec_seq, new_rec_seq):
            if seq == unique_seq:
                reconstructed_seq.append(new_seq)    
    return reconstructed_seq

def blast_analysis(df, db, pep=True, rec=False, generate_tmp_pep_df=False, generate_tmp_rec_df=False):
    """
    Perform BLAST analysis.
    This function runs a BLAST analysis on sequences from a DataFrame and generate output files used for mmseqs MSA search. 

    Parameters
    ----------
    df : pd.DataFrame
        Alphafold input DataFrame containing pdb id , receptor and peptides sequences.
    db : str
        Path to the BLAST database.
    pep : bool, optional
        Whether to perform analysis on peptide sequences (default True).
    rec : bool, optional
        Whether to perform analysis on nucleotide sequences (default False).
    generate_tmp_pep_df : bool, optional
        Whether to generate temporary output files for peptide sequences (default False).
    generate_tmp_rec_df : bool, optional
        Whether to generate temporary output files for nucleotide sequences (default False).

    Returns
    -------
    pd.DataFrame
        DataFrame containing additional information where for each pdb, the best sequence 
        with the highest coverage and length for the receptor and peptide.
    """
    print("Ruuning blast analysis")
    output_folder = './blast_analysis'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    blast_output = './blast_output'
    if not os.path.exists(blast_output):
        print(f"Blastp output folder is missing \n")


    query = []
    original_rec_seq = []
    original_pep_seq = []
    new_rec_seq = []
    new_pep_seq = []
    
    for full_rec, unique_rec, pep_seq, pdb in tqdm(zip(df["full_rec_seq"], df["unique_rec_seq"], df["pep_seq"], df["pdb"]), total=len(df)):
        tmp_rec_df = None  
        tmp_pep_df = None  
        new_rec=[]

        if pep:
            pep_csv_path = os.path.join(blast_output, f"out_pep_{pdb}.csv")
            tmp_pep_df = process_blast_results(pep_csv_path, pep_seq, db)
            if generate_tmp_pep_df:
                print(f"Generating blast analysis results of {pdb} peptide \n")
                output_file_pep = os.path.join(output_folder, f"tmp_pep_df_{pdb}.csv")
                tmp_pep_df.to_csv(output_file_pep, index=False)
            if tmp_pep_df is not None:
                best_pep = best_seq(tmp_pep_df)
                new_pep_seq.append(best_pep)
            else:
                new_pep_seq.append("") 

        if rec :
            for i, seq in enumerate(unique_rec):
                rec_csv_path = os.path.join(blast_output, f"out_rec_{pdb}_{i}.csv")
                tmp_rec_df = process_blast_results(rec_csv_path, seq, db)
                if generate_tmp_rec_df:
                    print(f"Generating blast analysis results of {pdb} receptor \n ")
                    output_file_rec = os.path.join(output_folder, f"tmp_rec_df_{pdb}_{i}.csv")
                    tmp_rec_df.to_csv(output_file_rec, index=False)
                
                if tmp_rec_df is None:
                    new_rec_seq.append("")
                else:
                    best_rec = best_seq(tmp_rec_df)
                    new_rec.append(best_rec)
            if full_rec is None :
                original_rec_seq.append(unique_rec)
                new_rec_seq.append(new_rec)
            else:
                original_rec_seq.append(full_rec)
                new_rec_seq.append(reconstruct(full_rec,unique_rec,new_rec))
        query.append(pdb)
        original_pep_seq.append(pep_seq)        
    
    
    print(f"Preparing mmseq csv file \n")
    data = {
        "pdb": query,
        "original_rec_seq": original_rec_seq,
        "original_pep_seq": original_pep_seq,
    }
    mmseq_df = pd.DataFrame(data)
    if pep :
        mmseq_df["new_pep_seq"]=new_pep_seq
    if rec:
        mmseq_df["new_rec_seq"]=new_rec_seq

    output_mmseq = os.path.join(output_folder, f"mmseq.csv")
    mmseq_df.to_csv(output_mmseq, index=False)
    
