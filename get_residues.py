import json
from load_matrix import find_indices, plot

protein_file = ('/home/kszyman/Downloads/fold_protein_glycan_ion_pdb_7bbv/'
                'fold_protein_glycan_ion_pdb_7bbv_full_data_0.json')
request_file = ('/home/kszyman/Downloads/fold_protein_glycan_ion_pdb_7bbv/'
                'fold_protein_glycan_ion_pdb_7bbv_job_request.json')


with open(request_file) as json_file:
    data = json.load(json_file)[0]

sequence = data['sequences'][0]['proteinChain']['sequence']
print(sequence)
print(len(sequence))
for i in range(100, 700, 100):
    probs, indices = find_indices(protein_file, i)
    plot(probs, indices)
outfile = open("corr.txt", 'w')
for ind in indices:
    print((ind[0], sequence[ind[0]]), (ind[1], sequence[ind[1]]), file=outfile)
