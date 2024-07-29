from pydca.plmdca import plmdca
from pydca.meanfield_dca import meanfield_dca
from pydca.sequence_backmapper import sequence_backmapper
from pydca.msa_trimmer import msa_trimmer
from pydca.contact_visualizer import contact_visualizer
from pydca.dca_utilities import dca_utilities
from Bio import SeqIO
from Bio import Seq
import os


msa_file = "/home/kszyman/Downloads/seqdump_padded.fasta"
refseq_file = '/home/kszyman/Downloads/sequence2.fasta'

# create MSATrimmer instance
trimmer = msa_trimmer.MSATrimmer(
    msa_file, biomolecule='protein',
    refseq_file=refseq_file,
)


"""input_file = '/home/kszyman/Downloads/seqdump.txt'
records = SeqIO.parse(input_file, 'fasta')
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment
output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
print(output_file)
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
"""
trimmed_data = trimmer.get_msa_trimmed_by_refseq(remove_all_gaps=False)

# write trimmed msa to file in FASTA format
trimmed_data_outfile = '/home/kszyman/Downloads/Trimmed.fa'
with open(trimmed_data_outfile, 'w') as fh:
    for seqid, seq in trimmed_data:
        fh.write('>{}\n{}\n'.format(seqid, seq))

plmdca_inst = plmdca.PlmDCA(
    trimmed_data_outfile,
    'rna',
    seqid = 0.8,
    lambda_h = 1.0,
    lambda_J = 20.0,
    num_threads = 10,
    max_iterations = 500,
)

# compute DCA scores summarized by Frobenius norm and average product corrected
plmdca_FN_APC = plmdca_inst.compute_sorted_FN_APC()

for site_pair, score in plmdca_FN_APC[:50]:
    print(site_pair, score)
