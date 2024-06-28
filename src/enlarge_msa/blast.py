import os 
import pandas as pd
import subprocess
import glob
import shutil
from tqdm import tqdm

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
    pdb_list = AF_seq["id"].tolist()
    
    for seq in AF_seq["sequence"].tolist():
        pep = seq.split(":")[-1]
        rec = ":".join(seq.split(":")[0:-1])  # Joining the receptor sequence parts back together
        pep_seq.append(pep)
        rec_seq.append(rec)
    
    if len(pep_seq) == len(rec_seq) == len(pdb_list):
        data = {
            "pdb": pdb_list,
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
    rec_file_path = os.path.join(folder_name, f"rec_{pdb}.fa")
    pep_file_path = os.path.join(folder_name, f"pep_{pdb}.fa")
    print(pdb)
    print(f"prepraring rec_{pdb} fasta file \n")
    print(f"prepraring pep_{pdb} fasta file \n")
    with open(rec_file_path, "w") as file:
        file.write(f">rec_{pdb}\n")
        file.write(str(rec_seq))
    with open(pep_file_path, "w") as file:
        file.write(f">pep_{pdb}\n")
        file.write(str(pep_seq))

def blasp_launch(pdb,fasta, db, output ,num_threads ,evalue):
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

    Returns
    -------
    None
    """
    
    BLASTP_BIN = shutil.which('blastp')
    if BLASTP_BIN is None:
        print('blastp could not be founded')

    command_line = [BLASTP_BIN,'-query',
                    fasta,'-db',db ,'-out',
                    output,'-outfmt', '10', '-num_threads', str(num_threads), '-evalue', str(evalue)]
    subprocess.Popen(command_line)

def blastp (df, db, evalue,num_threads , pep=True, rec=False):
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

    Returns
    -------
    None
    """

    for rec_seq, pep_seq, pdb in zip(df["rec_seq"], df["pep_seq"], df["pdb"]):
        blast_fasta_generator(rec_seq, pep_seq, pdb)
        if pep:
            fasta= f'pep_{pdb}.fa'
            output=f'out_pep_{pdb}.csv'
            blasp_launch(pdb, fasta, db, output ,str(num_threads),str(evalue))
        if rec: 
            fasta= f'rec_{pdb}.fa'
            output=f'out_rec_{pdb}.csv'
            blasp_launch(pdb, fasta, db, output ,str(num_threads),str(evalue))
        else : 
            print("No blast type selected")
