import subprocess
from Bio import SeqIO
from helper import reformat_result
import pandas as pd
    
def get_blast_result(inpname = 'uploaded.fasta', database = "AVP_BLASTDB", outname = "avp_blast_output.csv"):
    command = f'blastp|-db|assets/blastdb/{database}|-query|input/{inpname}|-out|output/{outname}|-max_target_seqs|1|-outfmt|6' 
    print(command.split('|'))
    subprocess.run(command.split('|'))

    blast_df = pd.read_csv(f'output/{outname}', sep = '\t', names ='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split())
    blast_df['original_index'] = blast_df['sseqid'].apply(lambda x: int(x.split('_')[0][3:]))
    original_data = pd.read_csv('assets/avp_data_v2.csv').reset_index()
    avp_merged_df = pd.merge(blast_df, original_data, left_on = 'original_index', right_on = 'index')
    reformat_result(avp_merged_df, rename_dict = {'qseqid':'Query sequence ID', 'sseqid': 'Matched sequence ID', 'pident': 'Percent identity', 'Class': 'Class', 'Source': 'Data source'}).to_csv(f'output/avp_blast_output.csv', index=False)

    anticov_result_df = pd.read_csv('output/anticov_result.csv')
    anticov_merged_df = avp_merged_df[avp_merged_df['qseqid'].isin(anticov_result_df['Sequence ID'])]
    reformat_result(anticov_merged_df, rename_dict = {'qseqid':'Query sequence ID', 'sseqid': 'Matched sequence ID', 'pident': 'Percent identity', 'Class': 'Class', 'Source': 'Data source'}).to_csv(f'output/anticov_blast_output.csv', index=False)