import pandas as pd 
import numpy as np
import copy
from tqdm.auto import tqdm 
from rdkit import Chem 
tqdm.pandas()


def fp_calculation(df , descriptor = "maccs"):
    from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect 
    from rdkit.Chem import MACCSkeys
    from rdkit.Chem import Descriptors
    moldf = copy.deepcopy(df)
    if (descriptor=='maccs'):
        maccs_fp = lambda mol: MACCSkeys.GenMACCSKeys(mol)  
        moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
        moldf['Descriptors'] = moldf.Mol.progress_apply(maccs_fp)
    elif (descriptor=='morgan'):  
        nBits = 2 * 1024
        radius= 2
        calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=False)
        moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
        moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)      
    elif (descriptor=='morgan4'):    
        nBits = 2 * 1024
        radius= 4
        calcfp = lambda mol: GetMorganFingerprintAsBitVect(mol,nBits=nBits,radius=radius, useFeatures=False, useChirality=False)
        moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
        moldf['Descriptors'] = moldf.Mol.progress_apply(calcfp)       
    else:
        print("The descriptor is not implemented")
    return moldf

def modi(df , outcome , descriptor = "maccs"):
    moldf = copy.deepcopy(df)
    import pandas as pd 
    import numpy as np
    from scipy.spatial import distance
    moldf = fp_calculation(moldf, descriptor)
    try:
        a = np.array(list(moldf['Descriptors']))
        dm = pd.DataFrame(distance.cdist(a, a, 'euclidean'))
        dm["nn_dist"] = dm.replace(0, np.nan).min()
        dm2 = moldf.merge(dm, left_index=True, right_index=True)    
        outcomeList = []
        for x in range(dm2.shape[0]):
            sel_index= dm2[x].index[dm2[x] == dm2["nn_dist"][x]].values
            outcomeList.append(dm2[outcome][dm2.index == sel_index[0]].to_list()[0])
        dm2["nn_outcome"] = outcomeList
        tot_agree = dm2[dm2[outcome] == dm2["nn_outcome"]].shape[0]
        tot = dm2[dm2[outcome] == dm2["nn_outcome"]].shape[1]
        modi = tot_agree/tot
    except:
        print("The descriptor is not implemented")
        modi = None
    return modi 