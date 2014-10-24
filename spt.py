from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np

suppl = Chem.SDMolSupplier('syk.sdf')

fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in suppl]

mat = []
for i in range(len(fps)):
    mat.append([DataStructs.FingerprintSimilarity(fps[i],fps[j]) for j in range(len(fps))])


similarity_matrix = np.array(mat)

edges = []
remain = set(range(1,len(fps)))

current = set([i for i,s in enumerate(similarity_matrix[0]) if s > 0.55 and i != 0])
processed = set([i for i,s in enumerate(similarity_matrix[0]) if s > 0.4 and i != 0]) - current

remain = remain - current - processed

for v in current:
    edges.append([0,v])


print current, processed
print remain


print "ID\tSmiles"
for i, mol in enumerate(suppl):
   print "{}\t{}".format(i, Chem.MolToSmiles(mol))

