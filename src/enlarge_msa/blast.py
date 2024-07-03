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

def process_csv(path):
    """
    Prepare a DataFrame from an Alphafold input CSV file.

    This function reads an Alphafold input CSV file to extract and prepare a DataFrame containing 
    separated peptide and receptor sequences.

    Parameters
    ----------
    path : str
        The path to the Alphafold input CSV file.

    Returns
    -------
    df : pandas.DataFrame or None
        DataFrame containing columns 'pdb', 'rec_seq', and 'pep_seq' if the CSV file is processed successfully.
        None if there is an issue with the CSV file.
    """
    AF_seq = pd.read_csv(path)
    
    pep_seq = []
    rec_seq = []
    original_rec_seq=[]
    pdb_list = AF_seq["id"].tolist()
    
    for seq in AF_seq["sequence"].tolist():
        pep = seq.split(":")[-1]
        rec = seq.split(":")[0:-1]
        if len(rec)!= len(set(rec)):
            original_rec_seq.append(rec)
        else : 
            original_rec_seq.append(None)
        pep_seq.append(pep)
        rec_seq.append(pd.unique(np.array(rec)).tolist())
    
    if len(pep_seq) == len(rec_seq) == len(pdb_list):
        data = {
            "pdb": pdb_list,
            "original_rec_seq":original_rec_seq,
            "rec_seq": rec_seq,
            "pep_seq": pep_seq,
        }
        
        df = pd.DataFrame(data)
        return df
    else:
        print("Error: problem with the csv.")
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
    folder_name = '.'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    pep_file_path = os.path.join(folder_name, f"pep_{pdb}.fa")
    print(f"Preparing pep_{pdb} FASTA file\n")
    
    with open(pep_file_path, "w") as file:
        file.write(f">pep_{pdb}\n")
        file.write(pep_seq)
        
    for i, seq in enumerate(rec_seq):
        rec_file_path = os.path.join(folder_name, f"rec_{pdb}_{i}.fa")
        print(f"Preparing rec_{pdb}_{i} FASTA file\n")
        
        with open(rec_file_path, "w") as file:
            file.write(f">rec_{pdb}_{i}\n")
            file.write(seq)
    


def blasp_launch(pdb, fasta, db, output, num_threads, evalue, semaphore):
    """
    Launch a BLASTP search.

    This function runs a BLASTP search using the specified parameters. It constructs the BLASTP command line
    and executes it as a subprocess.

    Parameters
    ----------
    pdb : str
        The PDB of the query associated with the input FASTA file.
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
            print('blastp could not be found')
            return

        command_line = [
            BLASTP_BIN, '-query', fasta, '-db', db, '-out', output,
            '-outfmt', '10', '-num_threads', str(num_threads), '-evalue', str(evalue)
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

    for rec_seq, pep_seq, pdb in zip(df["rec_seq"], df["pep_seq"], df["pdb"]):
        blast_fasta_generator(rec_seq, pep_seq, pdb)  
        if pep:
            tasks_pep.append((pdb, f'pep_{pdb}.fa', db, f'out_pep_{pdb}.csv', num_threads, evalue, semaphore))
        if rec:
            for i, seq in enumerate(rec_seq):
                tasks_rec.append((pdb, f'rec_{pdb}_{i}.fa', db, f'out_rec_{pdb}_{i}.csv', num_threads, evalue, semaphore))
                
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
        print('blastdbcmd could not be founded')
        
    command_line = [Blastdbcmd_BIN,'-db',db ,'-entry',
                    query]
    result = subprocess.run(command_line, capture_output=True, text=True)
    lines = result.stdout.splitlines()
    seq = lines[1] if len(lines) > 1 else ""

    return seq



def tmp_results(csv_path, seq, db):
    """
    Treat BLAST results and augment with coverage percentage and sequences.

    This function reads a CSV file of BLAST search results, calculates the coverage percentage 
    for each alignment, retrieves the corresponding sequences from the BLAST searched database, and adds 
    these sequences and their lengths to the DataFrame.

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
        A DataFrame containing the original BLAST search results, with additional columns for 
        coverage percentage, sequences, and sequence lengths.
    """
    
    tmp_df=None
    column_names = ["query", "subject_id", "% identity", "alignment_length", "mismatches", 
                "gap_open", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    coverage = []
    seqs = []
    tmp_df = pd.read_csv(csv_path, names=column_names)
    if tmp_df is not None :
        coverage = [(int(i) / len(seq)) * 100 for i in tmp_df["alignment_length"].tolist()]
        tmp_df["coverage %"] = coverage
        
        for query in tmp_df["subject_id"].tolist():
            new_seq = blastdbcmd(query, db)
            if new_seq:
                seqs.append(new_seq)
            else:
                seqs.append("")  # Handle case where no sequence is found
        
        tmp_df["sequence"] = seqs
        tmp_df["sequence_length"] = [len(i) for i in seqs]
    
    return tmp_df


def best_seq(tmp_df):
    """
    Identify and return the best sequence based on maximum coverage and sequence length.

    This function finds first the sequence with the highest coverage percentage in the provided DataFrame.
    If there are multiple sequences with the same highest coverage, it returns the longest sequence.

    Parameters
    ----------
    tmp_df : pd.DataFrame
        The DataFrame containing BLAST search results with additional columns for coverage percentage 
        and sequence information.

    Returns
    -------
    str
        The sequence with the highest coverage percentage and the longest length. Returns an empty 
        string if the DataFrame is empty.
    """
    
    if tmp_df.empty:
        return ""  # Return an empty string or handle appropriately

    max_coverage = tmp_df['coverage %'].max()
    max_coverage_df = tmp_df[tmp_df['coverage %'] == max_coverage]
    sequence_with_max_coverage_and_length = tmp_df.loc[max_coverage_df['sequence_length'].idxmax()]
    
    query = sequence_with_max_coverage_and_length["query"]
    best_id = sequence_with_max_coverage_and_length["subject_id"]
    new_seq = sequence_with_max_coverage_and_length["sequence"]

    return new_seq


def reconstruct(original_rec, rec_seq, new_rec):
    """
    Reconstruction of a sequence.

    This function helps in reconstructing an original sequence with new best retrieved ones 
    for the case of proteins with repetitive sequences. It replaces each segment in the 
    original sequence with new ones, in order to obtain the final full sequence.

    Parameters
    ----------
    original_rec : list of str
        The original sequence containing repetitive segments to be replaced.
    rec_seq : list of str
        A list of unique segments obtained from the original sequence.
    new_rec : list of str
        A list of new sequences that will replace the corresponding segments in the original one.

    Returns
    -------
    list of str
        The reconstructed sequence with segments replaced as specified.
        
    Examples
    --------
    >>> original_rec = ['A', 'B', 'A']
    >>> rec_seq = ['A', 'B']
    >>> new_rec = ['X', 'Y']
    >>> reconstruct(original_rec, rec_seq, new_rec)
    ['X', 'Y', 'X']
    """

    reconstructed_seq = []
    for seq in original_rec:
        for unique_seq, new_seq in zip(rec_seq, new_rec):
            if seq == unique_seq:
                reconstructed_seq.append(new_seq)    
    return reconstructed_seq

def blast_analysis(df, db, pep=True, rec=False, generate_tmp_pep_files=False, generate_tmp_rec_files=False):
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
    generate_tmp_pep_files : bool, optional
        Whether to generate temporary output files for peptide sequences (default False).
    generate_tmp_rec_files : bool, optional
        Whether to generate temporary output files for nucleotide sequences (default False).

    Returns
    -------
    pd.DataFrame
        DataFrame containing additional information where for each pdb, the best sequence 
        with the highest coverage and length for the receptor and peptide.
    """
    
    output_folder = 'blast_analysis'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    query = []
    original_rec_seq = []
    original_pep_seq = []
    new_rec_seq = []
    new_pep_seq = []
    
    for original_rec,rec_seq, pep_seq, pdb in tqdm(zip(df["original_rec_seq"],df["rec_seq"], df["pep_seq"], df["pdb"])):
        tmp_rec = None  
        tmp_pep = None  
        new_rec=[]

        if pep:
            csv_path = f"out_pep_{pdb}.csv"
            tmp_pep = tmp_results(csv_path, pep_seq, db)
            if generate_tmp_pep_files:
                output_file_pep = os.path.join(output_folder, f"tmp_pep_{pdb}.csv")
                tmp_pep.to_csv(output_file_pep, index=False)
            if tmp_pep is not None:
                best_pep = best_seq(tmp_pep)
                new_pep_seq.append(best_pep)
            else:
                new_pep_seq.append("") 

        #######     
        if rec :
            for i, seq in enumerate(rec_seq):
                csv_path = f"out_rec_{pdb}_{i}.csv"
                tmp_rec = tmp_results(csv_path, seq, db)
                if generate_tmp_rec_files:
                    output_file_rec = os.path.join(output_folder, f"tmp_rec_{pdb}_{i}.csv")
                    tmp_rec.to_csv(output_file_rec, index=False)
                
                if tmp_rec is None:
                    new_rec_seq.append("")
                else:
                    best_rec = best_seq(tmp_rec)
                    new_rec.append(best_rec)
                    
            if original_rec is None :
                original_rec_seq.append(rec_seq)
                new_rec_seq.append(new_rec)
            else:
                original_rec_seq.append(original_rec)
                new_rec_seq.append(reconstruct(original_rec,rec_seq,new_rec))
                
        query.append(pdb)
        original_pep_seq.append(pep_seq)        
    
    

    data = {
        "pdb": query,
        "original_rec_seq": original_rec_seq,
        "original_pep_seq": original_pep_seq,
       #"new_rec_seq": new_rec_seq,
       #"new_pep_seq": new_pep_seq
    }
    mmseq_df = pd.DataFrame(data)
    if pep :
        mmseq_df["new_pep_seq"]=new_pep_seq
    if rec:
        mmseq_df["new_rec_seq"]=new_rec_seq

    output_mmseq = os.path.join(output_folder, f"mmseq.csv")
    mmseq_df.to_csv(output_mmseq, index=False)
    




