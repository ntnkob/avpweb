import pycurl
import certifi
import subprocess
from io import BytesIO
import pandas as pd
from helper import reformat_result

buffer = BytesIO()

def structure_process(all_fasta):
    for seq_id, sequence in zip(all_fasta['id'], all_fasta['seq']):
        with open(f'output/predictedPDB/{seq_id}.pdb', 'w') as writehandle:
            writehandle.write(get_PDB_structure(sequence))

    command = f'foldseek easy-search output/predictedPDB assets/structuredb/avp_structdb output/structure_result.tsv tmp --max-seqs 1 --format-mode 4' 
    print(command.split(' '))
    subprocess.run(command.split(' '))

    structure_df = pd.read_csv(f'output/structure_result.tsv', sep = '\t')
    structure_df['original_index'] = structure_df['target'].apply(lambda x: int(x.split('_')[0][3:]))
    original_data = pd.read_csv('assets/avp_data_v2.csv').reset_index()

    avp_merged_df = pd.merge(structure_df, original_data, left_on = 'original_index', right_on = 'index')
    reformat_result(avp_merged_df, rename_dict = {'query':'Query sequence ID', 'target': 'Matched sequence ID', 'fident': 'Percent identity', 'Class': 'Class', 'Source': 'Data source'}).to_csv(f'output/avp_structure_result.csv', index=False)

    anticov_result_df = pd.read_csv('output/anticov_result.csv')
    anticov_merged_df = avp_merged_df[avp_merged_df['query'].isin(anticov_result_df['Sequence ID'])]
    reformat_result(avp_merged_df, rename_dict = {'query':'Query sequence ID', 'target': 'Matched sequence ID', 'fident': 'Percent identity', 'Subclass': 'Class', 'Source': 'Data source'}).to_csv(f'output/anticov_structure_result.csv', index=False)


def get_PDB_structure(sequence = 'IFKAIWSGIK'):
    c = pycurl.Curl()
    c.setopt(c.CAINFO, certifi.where())
    c.setopt(c.SSL_VERIFYPEER, 0)
    c.setopt(c.URL, 'https://api.esmatlas.com/foldSequence/v1/pdb/')
    c.setopt(c.POSTFIELDS, sequence)
    c.setopt(c.WRITEDATA, buffer)
    c.perform()
    c.close()

    body = buffer.getvalue()
    output = body.decode('utf-8')
    buffer.flush()
    buffer.seek(0)
    return output

print(get_PDB_structure())