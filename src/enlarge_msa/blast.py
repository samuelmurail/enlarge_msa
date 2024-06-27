import subprocess
import shutil

BLASTP_BIN = shutil.which('blastp')
if BLASTP_BIN is None:
    print('blastp could not be founded')
    BLAST_BIN = '/shared/software/conda/envs/blast-2.14.0/bin/blastp'


def blasp(seq, db, out_name):

    # Ecrire seq dans tmp.fa

    command_line = [BLASTP_BIN,'-query',
                    'tmp.fa','-out',
                    out_name,'-db',db]
    subprocess.call(command_line)