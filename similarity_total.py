#! bin/usr/env python

"""
Compute molecular similarities based on ap,morgan,cmorgan,fmorgan fingerprints
"""

import cPickle, gzip, numpy, copy, math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
#from matplotlib import cm
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB

#hard coded variables
bit_size = 2048 
radius = 1
output = open('strobi_similarity_total_2048','w')
output.write('reference compound ap morgan2 cmorgan2 fmorgan2\n')

# load reference compound and two test molecules
mols = []
compounds = []
for line in open('strobi_smile.txt'):
    mols.append(Chem.MolFromSmiles(line.split()[1]))
    compounds.append(line.split()[0])
def get_similarity(): # get similarities on the first molecule in compound group
    # precalculate fingerprints for reference compound
    ref_morgan2 = AllChem.GetMorganFingerprintAsBitVect(mols[0],radius,bit_size)
    ref_cmorgan2 = AllChem.GetMorganFingerprint(mols[0],radius)
    ref_fmorgan2 = AllChem.GetMorganFingerprintAsBitVect(mols[0], radius,bit_size, useFeatures = True)
    ref_ap = Pairs.GetAtomPairFingerprint(mols[0])
    # precalculate fingerprints and bit information for test molecules
    total_sims = ''
    fps_morgan2 = []
    fps_cmorgan2 = []
    fps_fmorgan2 = []
    fps_ap = []
    info_morgan2 = []
    info_cmorgan2 = []
    info_fmorgan2 = []
    num_mols = len(mols) - 1
    reference = compounds[0]
    del compounds[0] 
    del mols[0] #remove reference cmp from list
    for m in mols:
        info = {}
        fps_morgan2.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, bit_size,  bitInfo = info))
        info_morgan2.append(info)
        info = {}
        fps_cmorgan2.append(AllChem.GetMorganFingerprint(m, radius, bitInfo=info))
        info_cmorgan2.append(info)
        info = {}
        fps_fmorgan2.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, bit_size, useFeatures=True, bitInfo=info))
        info_fmorgan2.append(info)
        fps_ap.append(Pairs.GetAtomPairFingerprint(m))
    ## calculate similarities
    for i,m in enumerate(mols):
        ap_simil = DataStructs.DiceSimilarity(ref_ap, fps_ap[i])
        morgan2_simil = DataStructs.DiceSimilarity(ref_morgan2, fps_morgan2[i])
        cmorgan2_simil = DataStructs.DiceSimilarity(ref_cmorgan2, fps_cmorgan2[i])
        fmorgan2_simil = DataStructs.DiceSimilarity(ref_fmorgan2, fps_fmorgan2[i])
        sims =str(reference)+' '+ str(compounds[i].rstrip())+' '+ str(ap_simil)+' '+str(morgan2_simil)+' '+str(cmorgan2_simil)+' '+str(fmorgan2_simil)+'\n'
        total_sims += sims
    return total_sims
for i in range(len(compounds)):
    output.write(str(get_similarity()))   
output.close()
