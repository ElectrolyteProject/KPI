import os
import pandas as pd
import numpy as np
import csv
import math
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import Fragments
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

# Molecular Weight
def MolWt(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.MolWt(mol)

# Heavy Atom Count
def Heavy_atoms(mol):
    return Lipinski.HeavyAtomCount(mol)

# Carbon Atom Count
def count_C_atom(smiles):
    if "Cl" in smiles or "cl" in smiles:
        return smiles.count('C') + smiles.count('c') - count_Cl_atom(smiles)
    else:
        return smiles.count('C') + smiles.count('c')

# Oxygen Atom Count
def count_O_atom(smiles):
    return smiles.count("O") + smiles.count("o")

# Oxygen to Carbon Atom Ratio
def count_O_C_atom(smiles):
    return count_O_atom(smiles) / (count_C_atom(smiles) + 0.001)

# Heteroatom Count
def count_Het_atom(mol):
    return Lipinski.NumHeteroatoms(mol)

# Heteroatom to Carbon Atom Ratio
def count_Het_C_atom(mol, smiles):
    return count_Het_atom(mol) / (count_C_atom(smiles) + 0.001)

# Double Bond Count
def double_bond(smiles):
    return smiles.count('=')

# Triple Bond Count
def triple_bond(smiles):
    return smiles.count('#')

# Hydrogen Bond Donors
def Donor(mol):
    return Lipinski.NumHDonors(mol)

# Hydrogen Bond Acceptors
def Accept(mol):
    return Lipinski.NumHAcceptors(mol)

# Rotatable Bond Count
def Rot(mol):
    return Lipinski.NumRotatableBonds(mol)

# Ring Count
def Rings(smiles):
    mol = Chem.MolFromSmiles(smiles)
    ri = mol.GetRingInfo()
    return ri.NumRings()

# Largest Ring Atom Count
def Nring(molecule):
    rings = molecule.GetRingInfo().AtomRings()
    if not rings:  # Check if there are any rings
        return 0  # If no rings, return 0
    largest_ring = max(rings, key=len)  # Find the largest ring
    return len(largest_ring)  # Return the number of atoms in the largest ring

# Aliphatic Carbocyclic Ring Count
def AlCR(mol):
    return Lipinski.NumAliphaticCarbocycles(mol)

# Aliphatic Heterocyclic Ring Count
def AlHR(mol):
    return Lipinski.NumAliphaticHeterocycles(mol)

# Aliphatic Ring Count
def AlR(mol):
    return Lipinski.NumAliphaticRings(mol)

# Aromatic Carbocyclic Ring Count
def ArCR(mol):
    return Lipinski.NumAromaticCarbocycles(mol)

# Aromatic Heterocyclic Ring Count
def ArHR(mol):
    return Lipinski.NumAromaticHeterocycles(mol)

# Aromatic Ring Count
def ArR(mol):
    return Lipinski.NumAromaticRings(mol)

# Saturated Carbocyclic Ring Count
def SCR(mol):
    return Lipinski.NumSaturatedCarbocycles(mol)

# Saturated Heterocyclic Ring Count
def SHR(mol):
    return Lipinski.NumSaturatedHeterocycles(mol)

# Saturated Ring Count
def SR(mol):
    return Lipinski.NumSaturatedRings(mol)

def find_largest_ring(molecule):
    rings = molecule.GetRingInfo().AtomRings()
    if not rings:
        return []
    largest_ring = max(rings, key=len)
    return largest_ring

def find_longest_carbon_chain(molecule):
    G = rdmolops.GetAdjacencyMatrix(molecule)
    graph = nx.Graph(G)
    carbon_atoms = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() == 'C']
    longest_chain = []
    for start_atom in carbon_atoms:
        for end_atom in carbon_atoms:
            for path in nx.all_simple_paths(graph, start_atom, end_atom):
                if set(path).issubset(carbon_atoms) and len(path) > len(longest_chain):
                    longest_chain = path
    return longest_chain

def count_carbon_side_chains(molecule, main_chain):
    side_chain_count = 0
    for atom_idx in main_chain:
        atom = molecule.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in main_chain:
                side_chain_count += 1
    return side_chain_count

# Carbon Side Chain Count
def Bran(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    largest_ring = find_largest_ring(molecule)
    main_chain = largest_ring if largest_ring else find_longest_carbon_chain(molecule)
    side_chains = count_carbon_side_chains(molecule, main_chain)
    return side_chains

# Hydroxyl Group Count
def OH(mol):
    return Fragments.fr_Al_OH(mol) + Fragments.fr_Ar_OH(mol)

# Aliphatic Hydroxyl Group Count
def Al_OH(mol):
    return Fragments.fr_Al_OH(mol)

# Aromatic Hydroxyl Group Count
def Ar_OH(mol):
    return Fragments.fr_Ar_OH(mol)

# Carboxyl Group Count
def COO(mol):
    return Fragments.fr_COO(mol)

# Aliphatic Carboxyl Group Count
def Al_COO(mol):
    return Fragments.fr_Al_COO(mol)

# Aromatic Carboxyl Group Count
def Ar_COO(mol):
    return Fragments.fr_Ar_COO(mol)

# Carbonyl Oxygen Count
def C_O(mol):
    return Fragments.fr_C_O(mol)

# Carbonyl Oxygen Count (excluding Carboxyl)
def C_O_noCOO(mol):
    return Fragments.fr_C_O_noCOO(mol)

# Ether Count
def ether(mol):
    return Fragments.fr_ether(mol)

# Ester Count
def ester(mol):
    return Fragments.fr_ester(mol)

# Halogen Count
def halogen(mol):
    return Fragments.fr_halogen(mol)

# Benzene Ring Count
def benzene(mol):
    return Fragments.fr_benzene(mol)

# Aromatic Nitrogen Count
def Ar_N(mol):
    return Fragments.fr_Ar_N(mol)

# Nitrogen Atoms in Aromatic Systems
def ArN(mol):
    return Fragments.fr_ArN(mol)

# Aromatic Amine Count
def Ar_NH(mol):
    return Fragments.fr_Ar_NH(mol)

# Imine Count
def Imine(mol):
    return Fragments.fr_Imine(mol)

# Primary Amine Count
def NH2(mol):
    return Fragments.fr_NH2(mol)

# Secondary Amine Count
def NH1(mol):
    return Fragments.fr_NH1(mol)

# Tertiary Amine Count
def NH0(mol):
    return Fragments.fr_NH0(mol)

# Thiol Count
def SH(mol):
    return Fragments.fr_SH(mol)

# Aldehyde Count
def aldehyde(mol):
    return Fragments.fr_aldehyde(mol)

# Amide Count
def amide(mol):
    return Fragments.fr_amide(mol)

# Aniline Count
def aniline(mol):
    return Fragments.fr_aniline(mol)

# Phenol Count
def phenol(mol):
    return Fragments.fr_phenol(mol)

# Epoxide Count
def epoxide(mol):
    return Fragments.fr_epoxide(mol)

# Furan Ring Count
def furan(mol):
    return Fragments.fr_furan(mol)

# Piperidine Ring Count
def piperdine(mol):
    return Fragments.fr_piperdine(mol)

# Pyridine Ring Count
def pyridine(mol):
    return Fragments.fr_pyridine(mol)

# Lactone Ring Count
def lactone(mol):
    return Fragments.fr_lactone(mol)

# Ketone Count
def ketone(mol):
    return Fragments.fr_ketone(mol)

# Nitrile Count
def nitrile(mol):
    return Fragments.fr_nitrile(mol)

# Sulfone Count
def sulfone(mol):
    return Fragments.fr_sulfone(mol)

# Hydrogen Atom Count
def count_H_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    supp_H = AllChem.AddHs(mol)
    smiles_with_H = Chem.MolToSmiles(supp_H)
    return smiles_with_H.count("H") + smiles_with_H.count("h")

# Nitrogen Atom Count
def count_N_atom(smiles):
    return smiles.count('N') + smiles.count('n')

# Fluorine Atom Count
def count_F_atom(smiles):
    return smiles.count("F") + smiles.count("f")

# Silicon Atom Count
def count_Si_atom(smiles):
    return smiles.count('Si') + smiles.count('si')

# Phosphorus Atom Count
def count_P_atom(smiles):
    return smiles.count("P") + smiles.count("p")

# Sulfur Atom Count
def count_S_atom(smiles):
    if "Si" in smiles or "si" in smiles:
        return smiles.count('S') + smiles.count('s') - count_Si_atom(smiles)
    else:
        return smiles.count('S') + smiles.count('s')

# Chlorine Atom Count
def count_Cl_atom(smiles):
    return smiles.count("Cl") + smiles.count("cl")

# Bromine Atom Count
def count_Br_atom(smiles):
    return smiles.count('Br') + smiles.count('br')

# Iodine Atom Count
def count_I_atom(smiles):
    if "Si" in smiles or "si" in smiles:
        return smiles.count("I") + smiles.count("i") - count_Si_atom(smiles)
    else:
        return smiles.count("I") + smiles.count("i")
    
# Total Atom Count
def count_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    supp_H = AllChem.AddHs(mol)
    smiles = Chem.MolToSmiles(supp_H)
    cou = 0
    for i in smiles:
        if i.isalpha():
            cou += 1
    return cou

# Average Electronegativity
def avg_pauling_electronegativity(smiles):
    pelec_H = 2.20
    pelec_C = 2.55
    pelec_N = 3.04
    pelec_O = 3.44
    pelec_F = 3.98
    pelec_Si = 1.90
    pelec_P = 2.19
    pelec_S = 2.58
    pelec_Cl = 3.16
    pelec_Br = 2.96
    pelec_I = 2.66
    return (pelec_H * count_H_atom(smiles) +
            pelec_C * count_C_atom(smiles) + 
            pelec_N * count_N_atom(smiles) +
            pelec_O * count_O_atom(smiles) +
            pelec_F * count_F_atom(smiles) + 
            pelec_Si * count_Si_atom(smiles) +
            pelec_P * count_P_atom(smiles) +
            pelec_S * count_S_atom(smiles) + 
            pelec_Cl * count_Cl_atom(smiles) +
            pelec_Br * count_Br_atom(smiles) +
            pelec_I * count_I_atom(smiles)
            ) / count_atom(smiles)

# Average First Ionization Energy
def avg_ionization_energies(smiles):
    ioni_H = 13.598433
    ioni_C = 11.26030
    ioni_N = 14.5341
    ioni_O = 13.61805
    ioni_F = 17.4228
    ioni_Si = 8.15168
    ioni_P = 10.48669
    ioni_S = 10.36001
    ioni_Cl = 12.96763
    ioni_Br = 11.8138
    ioni_I = 10.45126
    return (ioni_H * count_H_atom(smiles) +
            ioni_C * count_C_atom(smiles) + 
            ioni_N * count_N_atom(smiles) +
            ioni_O * count_O_atom(smiles) +
            ioni_F * count_F_atom(smiles) + 
            ioni_Si * count_Si_atom(smiles) +
            ioni_P * count_P_atom(smiles) +
            ioni_S * count_S_atom(smiles) + 
            ioni_Cl * count_Cl_atom(smiles) +
            ioni_Br * count_Br_atom(smiles) +
            ioni_I * count_I_atom(smiles)
            ) / count_atom(smiles)

# Average Electron Affinity
def avg_electron_affinity(smiles):
    ea_H = 0.75420375
    ea_C = 1.262118
    ea_N = -0.07
    ea_O = 1.4611120
    ea_F = 3.4011895
    ea_Si = 1.3895220
    ea_P = 0.7465
    ea_S = 2.0771029
    ea_Cl = 3.612724
    ea_Br = 3.3635880
    ea_I = 3.059038
    return (ea_H * count_H_atom(smiles) +
            ea_C * count_C_atom(smiles) + 
            ea_N * count_N_atom(smiles) +
            ea_O * count_O_atom(smiles) +
            ea_F * count_F_atom(smiles) + 
            ea_Si * count_Si_atom(smiles) +
            ea_P * count_P_atom(smiles) +
            ea_S * count_S_atom(smiles) + 
            ea_Cl * count_Cl_atom(smiles) +
            ea_Br * count_Br_atom(smiles) +
            ea_I * count_I_atom(smiles)
            ) / count_atom(smiles)
    
# Maximum Absolute Partial Charge
def MaxAPC(mol):
    return Descriptors.MaxAbsPartialCharge(mol)

# Minimum Absolute Partial Charge
def MinAPC(mol):
    return Descriptors.MinAbsPartialCharge(mol)

# Maximum Partial Charge
def MaxPC(mol):
    return Descriptors.MaxPartialCharge(mol)

# Minimum Partial Charge
def MinPC(mol):
    return Descriptors.MinPartialCharge(mol)

# Valence Electron Count
def ValE(mol):
    return Descriptors.NumValenceElectrons(mol)

def extract_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    raw_features = [
        MolWt(smiles), 
        Heavy_atoms(mol), 
        count_C_atom(smiles),
        count_O_atom(smiles),
        count_O_C_atom(smiles),
        count_Het_atom(mol),
        count_Het_C_atom(mol, smiles),
        
        double_bond(smiles),
        triple_bond(smiles),
        Donor(mol),
        Accept(mol),
        Rot(mol),
        Rings(smiles),
        Nring(mol),
        AlCR(mol),
        AlHR(mol),
        AlR(mol),
        ArCR(mol),
        ArHR(mol),
        ArR(mol),
        SCR(mol),
        SHR(mol),
        SR(mol),
        Bran(smiles),
        
        OH(mol),
        Al_OH(mol),
        Ar_OH(mol),
        COO(mol),
        Al_COO(mol),
        Ar_COO(mol),
        C_O(mol),
        C_O_noCOO(mol),
        ether(mol),
        ester(mol),
        halogen(mol),
        benzene(mol),
        Ar_N(mol),
        ArN(mol),
        Ar_NH(mol),
        Imine(mol),
        NH2(mol),
        NH1(mol),
        NH0(mol),
        SH(mol),
        aldehyde(mol),
        amide(mol),
        aniline(mol),
        phenol(mol),
        epoxide(mol),
        furan(mol),
        piperdine(mol),
        pyridine(mol),
        lactone(mol),
        ketone(mol),
        nitrile(mol),
        sulfone(mol),
        
        avg_pauling_electronegativity(smiles),
        avg_ionization_energies(smiles),
        avg_electron_affinity(smiles),
        MaxAPC(mol),
        MinAPC(mol),
        MaxPC(mol),
        MinPC(mol),
        ValE(mol),
        ]
    
    features = [0 if math.isnan(x) or math.isinf(x) else x for x in raw_features]
    return features

def feature_names():
    function_names = [
    "Molwt", 
    "#Heavy", 
    "#C",
    "#O",
    "#O/#C",
    "#Het",
    "#Het/#C",
    "#R=R",
    "#R#R",
    "#Donor",
    "#Accept",
    "#Rot",
    "#Ring",
    "#Nring",
    "#AlCR",
    "#AlHR",
    "#AlR",
    "#ArCR",
    "#ArHR",
    "#ArR",
    "#SCR",
    "#SHR",
    "#SR",
    "#Bran",
    "#OH",
    "#AlOH",
    "#ArOH",
    "#COO",
    "#AlCOO",
    "#ArCOO",
    "#C=O",
    "#C=O\COO",
    "#Ether",
    "#Ester",
    "#Halogen",
    "#Benzene",
    "#Ar-N",
    "#ArN",
    "#ArNH",
    "#Imine",
    "#NH2",
    "#NH1",
    "#NH0",
    "#SH",
    "#Aldehyde",
    "#Amide",
    "#Aniline",
    "#Phenol",
    "#Epoxide",
    "#Furan",
    "#Piperdine",
    "#Pyridine",
    "#Lactone",
    "#Ketone",
    "#Nitrile",
    "#Sulfone",
    "AvgX",
    "AvgI",
    "AvgA",
    "MaxAPC",
    "MinAPC",
    "MaxPC",
    "MinPC",
    "ValE"
    ]
    return function_names
def plot_pearson(data, cmap, folder_name):
    """
    Plot Pearson correlation heatmap and save feature correlations to a CSV file.

    Parameters:
    data (DataFrame): Input data containing features.
    cmap (str): Colormap for the heatmap.
    folder_name (str): Folder name to save the CSV file.
    """
    # Calculate Pearson correlation matrix
    corr_matrix = data.corr(method='pearson')
    
    # Plot heatmap
    sns.heatmap(corr_matrix, annot=False, cmap=cmap, vmin=-1, vmax=1)
    plt.show()
    
    # Extract correlations with the last column and sort by absolute value
    corr_last_column = corr_matrix.iloc[-1, :-1]
    corr_sorted = corr_last_column.abs().sort_values(ascending=False)
    
    # Create folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    # Save sorted correlations to CSV
    with open(f"./{folder_name}/1_heat_map_features_screening.csv", 'w') as f:
        writer = csv.writer(f)
        corr_sorted_values = [round(num, 3) for num in corr_sorted.values]
        writer.writerow(corr_sorted.index)
        writer.writerow(corr_sorted_values)

def save_data_to_csv(input_path, output_path, props):
    """
    Save molecular data and features to a CSV file.

    Parameters:
    input_path (str): Path to the input CSV file.
    output_path (str): Path to save the output CSV file.
    props (str): Property column name.
    """
    data = pd.read_csv(input_path)

    SMILES = data.loc[data["SMILES"].notnull(), "SMILES"].values
    props_list = data.loc[data[props].notnull(), props].values

    X = []
    y = []
    
    feature_name = feature_names()
    for smi, pr in zip(SMILES, props_list):
        try:
            X.append(extract_features(smi))
            y.append(pr)
        except Exception as e:
            print(f"Error processing SMILES {smi}: {e}")
    
    # Ensure data lengths match
    if len(SMILES) != len(X) or len(SMILES) != len(y):
        raise ValueError(f"SMILES: {len(SMILES)}, X: {len(X)}, and y: {len(y)} must have the same length")
    
    # Check feature name count matches feature count in X
    if len(feature_name) != len(X[0]):
        raise ValueError("The number of feature names must match the number of features in X")

    # Create DataFrame
    df = pd.DataFrame(X, columns=feature_name)
    df.insert(0, props, y)
    df.insert(0, 'SMILES', SMILES)

    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"The output file has {df.shape[0]} Molecules (rows) and {df.shape[1]} Features (columns).")

def save_data_to_csv_knoledge(input_path, output_path, props1, props2):
    """
    Save molecular data, features, and knowledge to a CSV file.

    Parameters:
    input_path (str): Path to the input CSV file.
    output_path (str): Path to save the output CSV file.
    props1 (str): First property column name.
    props2 (str): Second property column name.
    """
    data = pd.read_csv(input_path)

    SMILES = data.loc[data["SMILES"].notnull(), "SMILES"].values
    props_list = data.loc[data[props1].notnull(), props1].values

    X = []
    y = []
    
    feature_name = feature_names()
    for smi, pr in zip(SMILES, props_list):
        try:
            X.append(extract_features(smi))
            y.append(pr)
        except Exception as e:
            print(f"Error processing SMILES {smi}: {e}")
    
    # Ensure data lengths match
    if len(SMILES) != len(X) or len(SMILES) != len(y):
        raise ValueError(f"SMILES: {len(SMILES)}, X: {len(X)}, and y: {len(y)} must have the same length")
    
    # Check feature name count matches feature count in X
    if len(feature_name) != len(X[0]):
        raise ValueError("The number of feature names must match the number of features in X")

    # Create DataFrame
    df = pd.DataFrame(X, columns=feature_name)
    df.insert(0, props2, y)
    df.insert(0, 'SMILES', SMILES)

    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"The output file has {df.shape[0]} Molecules (rows) and {df.shape[1]} Features (columns).")

def create_folder(folder_name):
    """
    Check if a folder with the specified name exists. If not, create it.

    Parameters:
    folder_name (str): The name of the folder to check and create if it doesn't exist.
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Folder '{folder_name}' created.")
    else:
        print(f"Folder '{folder_name}' already exists.")
        