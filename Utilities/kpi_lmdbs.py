import os
import pickle
import lmdb
import pandas as pd
import numpy as np
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import RDLogger
from multiprocessing import Pool

# Disable RDKit logging
RDLogger.DisableLog('rdApp.*')

# Suppress warnings
import warnings
warnings.filterwarnings(action='ignore')

def smi2_2Dcoords(smi):
    """
    Convert SMILES to 2D coordinates.

    Parameters:
    smi (str): SMILES string.

    Returns:
    np.ndarray: 2D coordinates of the molecule.
    """
    mol = Chem.MolFromSmiles(smi)
    mol = AllChem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    coordinates = mol.GetConformer().GetPositions().astype(np.float32)
    assert len(mol.GetAtoms()) == len(coordinates), f"2D coordinates shape is not aligned with {smi}"
    return coordinates

def smi2_3Dcoords(smi, cnt):
    """
    Convert SMILES to 3D coordinates with multiple conformations.

    Parameters:
    smi (str): SMILES string.
    cnt (int): Number of conformers to generate.

    Returns:
    list: List of 3D coordinates of the molecule.
    """
    mol = Chem.MolFromSmiles(smi)
    mol = AllChem.AddHs(mol)
    coordinate_list = []

    for seed in range(cnt):
        try:
            res = AllChem.EmbedMolecule(mol, randomSeed=seed)
            if res == 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                    coordinates = mol.GetConformer().GetPositions()
                except:
                    print("Failed to generate 3D, using 2D instead")
                    coordinates = smi2_2Dcoords(smi)
            elif res == -1:
                mol_tmp = Chem.MolFromSmiles(smi)
                AllChem.EmbedMolecule(mol_tmp, maxAttempts=5000, randomSeed=seed)
                mol_tmp = AllChem.AddHs(mol_tmp, addCoords=True)
                try:
                    AllChem.MMFFOptimizeMolecule(mol_tmp)
                    coordinates = mol_tmp.GetConformer().GetPositions()
                except:
                    print("Failed to generate 3D, using 2D instead")
                    coordinates = smi2_2Dcoords(smi)
        except:
            print("Failed to generate 3D, using 2D instead")
            coordinates = smi2_2Dcoords(smi)

        assert len(mol.GetAtoms()) == len(coordinates), f"3D coordinates shape is not aligned with {smi}"
        coordinate_list.append(coordinates.astype(np.float32))

    return coordinate_list

def inner_smi2coords(content):
    """
    Convert SMILES to 2D and 3D coordinates.

    Parameters:
    content (tuple): Tuple containing SMILES, target, and additional knowledge.

    Returns:
    bytes: Serialized dictionary containing molecular data.
    """
    smi = content[0]
    target = content[1]
    knowledge = content[2:]
    cnt = 10  # Number of conformers: 10 3D + 1 2D

    mol = Chem.MolFromSmiles(smi)
    if len(mol.GetAtoms()) > 400:
        coordinate_list = [smi2_2Dcoords(smi)] * (cnt + 1)
        print(f"Atom count > 400, using 2D coordinates for {smi}")
    else:
        coordinate_list = smi2_3Dcoords(smi, cnt)
        coordinate_list.append(smi2_2Dcoords(smi).astype(np.float32))

    mol = AllChem.AddHs(mol)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return pickle.dumps({'atoms': atoms, 'coordinates': coordinate_list, 'mol': mol, 'smi': smi, 'target': target, 'knowledge': knowledge}, protocol=-1)

def smi2coords(content):
    """
    Wrapper function to handle SMILES conversion and error handling.

    Parameters:
    content (tuple): Tuple containing SMILES, target, and additional knowledge.

    Returns:
    bytes or None: Serialized dictionary containing molecular data or None if conversion fails.
    """
    try:
        return inner_smi2coords(content)
    except Exception as e:
        print(f"Failed to convert SMILES: {content[0]}, Error: {e}")
        return None

def write_lmdb(data_set, inpath='./', outpath='./', nthreads=16):
    """
    Write molecular data to an LMDB database.

    Parameters:
    data_set (str): Dataset name.
    inpath (str): Input CSV file path.
    outpath (str): Output directory path.
    nthreads (int): Number of threads for parallel processing.
    """
    df = pd.read_csv(os.path.join(inpath))
    content_list = list(zip(*[df[c].values.tolist() for c in df]))

    os.makedirs(outpath, exist_ok=True)
    output_name = os.path.join(outpath, f'{data_set}.lmdb')

    if os.path.exists(output_name):
        os.remove(output_name)

    env = lmdb.open(
        output_name,
        subdir=False,
        readonly=False,
        lock=False,
        readahead=False,
        meminit=False,
        max_readers=1,
        map_size=int(100e9),
    )
    txn = env.begin(write=True)
    
    with Pool(nthreads) as pool:
        for i, result in enumerate(tqdm(pool.imap(smi2coords, content_list))):
            if result is not None:
                txn.put(f'{i}'.encode("ascii"), result)
    
    txn.commit()
    env.close()

def knowledge2lmdbs(data_set, inpath, outpath):
    """
    Convert knowledge data to LMDB and calculate statistics.

    Parameters:
    data_set (str): Dataset name.
    inpath (str): Input CSV file path.
    outpath (str): Output directory path.
    """
    write_lmdb(data_set, inpath, outpath, nthreads=8)

    df = pd.read_csv(inpath)
    mean_value = df.iloc[:, 1].mean()
    std_value = df.iloc[:, 1].std()

    print(f'{data_set} mean: {mean_value}')
    print(f'{data_set} std: {std_value}')
