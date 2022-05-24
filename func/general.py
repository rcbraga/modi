import pandas as pd 
import numpy as np
import copy
from tqdm.auto import tqdm 
from rdkit import Chem 
tqdm.pandas()


def fp_calculation(df , descriptors = "maccs"):
    moldf = copy.deepcopy(df)
    from rdkit.Chem import MACCSkeys
    from rdkit.Chem import Descriptors
    maccs_fp = lambda mol: MACCSkeys.GenMACCSKeys(mol)  
    moldf['Mol'] = moldf.SMILES.apply(Chem.MolFromSmiles)
    moldf['Descriptors'] = moldf.Mol.progress_apply(maccs_fp)
    return moldf

def modi(df , outcome , descriptors = "maccs"):
    moldf = copy.deepcopy(df)
    import pandas as pd 
    import numpy as np
    from scipy.spatial import distance
    moldf = fp_calculation(moldf, descriptors)
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
    return tot_agree/tot  