import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
import logging

# Configure logging
log_file_path = 'process.log'
logging.basicConfig(level=logging.INFO, filename=log_file_path, filemode='w', format='%(message)s')

def extract_smiles_from_xyz(file_path):
    """
    Extracts the SMILES string from an XYZ file.

    Parameters:
    - file_path: str, path to the XYZ file.

    Returns:
    - str, the extracted SMILES string.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        smiles_line = lines[-2].strip().split()[0]  # Get the second last line and extract the first element
        return smiles_line

def convert_xyz_to_smiles(xyz_folder, output_file):
    """
    Converts all XYZ files in a folder to SMILES strings and saves them to an output file.

    Parameters:
    - xyz_folder: str, path to the folder containing XYZ files.
    - output_file: str, path to the output file to save SMILES strings.
    """
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(xyz_folder):
            if filename.endswith('.xyz'):
                file_path = os.path.join(xyz_folder, filename)
                smiles = extract_smiles_from_xyz(file_path)
                outfile.write(f"{smiles}\n")
    print(f"SMILES extraction complete. Results stored in {output_file}.")

def read_smiles_from_csv(file_path):
    """
    Reads SMILES and associated data from a CSV file.

    Parameters:
    - file_path: str, path to the CSV file.

    Returns:
    - DataFrame: DataFrame containing the data from the CSV file.
    """
    return pd.read_csv(file_path)

def has_active_hydrogen(smiles):
    """
    Checks if a molecule contains active hydrogen (hydrogen in hydroxyl, carboxyl, amine groups).

    Parameters:
    - smiles (str): SMILES string of the molecule.

    Returns:
    - bool: True if the molecule contains active hydrogen, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    # Check for hydroxyl (R-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return True

    # Check for primary amine (R-NH2), secondary amine (R-NHR'), tertiary amine (R-NR'R'')
    amine_patterns = [Chem.MolFromSmarts('[NX3;H2]'), Chem.MolFromSmarts('[NX3;H1]'), Chem.MolFromSmarts('[NX3;H0]')]
    for pattern in amine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True

    return False

def filter_initial_molecules(df):
    """
    Filters molecules based on allowed elements, heavy atoms, molecular weight, and active hydrogen.

    Parameters:
    - df: DataFrame, DataFrame containing the molecules to be filtered.

    Returns:
    - DataFrame: Filtered DataFrame.
    """
    allowed_elements = set(['H', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'])

    initial_count = len(df)
    logging.info(f"Initial number of molecules: {initial_count}")

    filtered_df = pd.DataFrame()

    for index, row in df.iterrows():
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            if elements.issubset(allowed_elements) and mol.GetNumHeavyAtoms() <= 30 and rdMolDescriptors.CalcExactMolWt(mol) <= 600 and not has_active_hydrogen(smiles):
                filtered_df = filtered_df.append(row, ignore_index=True)
    
    after_initial_filter_count = len(filtered_df)
    logging.info(f"After filtering active hydrogens, molecular weight <= 600, and heavy atoms <= 30: {after_initial_filter_count} molecules")

    return filtered_df

def filter_molecules(df, mp_threshold, bp_threshold, fp_threshold):
    """
    Further filters molecules based on melting point, boiling point, and flash point, and prints intermediate results.

    Parameters:
    - df: DataFrame, DataFrame containing the molecules to be filtered.
    - mp_threshold: float, melting point threshold.
    - bp_threshold: float, boiling point threshold.
    - fp_threshold: float, flash point threshold.

    Returns:
    - DataFrame: Filtered DataFrame.
    """
    # Filter based on melting point (mp)
    mp_filtered_df = df[df['mps_pred'] < mp_threshold]
    mp_filtered_count = len(mp_filtered_df)
    print(f"After filtering by mps_pred < {mp_threshold}K ({mp_threshold-273.15}°C): {mp_filtered_count} molecules")
    logging.info(f"After filtering by mps_pred < {mp_threshold}K: {mp_filtered_count} molecules")

    # Filter based on boiling point (bp)
    bp_filtered_df = mp_filtered_df[mp_filtered_df['bps_pred'] > bp_threshold]
    bp_filtered_count = len(bp_filtered_df)
    print(f"After filtering by bps_pred > {bp_threshold}K ({bp_threshold-273.15}°C): {bp_filtered_count} molecules")
    logging.info(f"After filtering by bps_pred > {bp_threshold}K: {bp_filtered_count} molecules")

    # Filter based on flash point (fp)
    final_filtered_df = bp_filtered_df[bp_filtered_df['fps_pred'] > fp_threshold]
    final_filtered_count = len(final_filtered_df)
    print(f"After filtering by fps_pred > {fp_threshold}K ({fp_threshold-273.15}°C): {final_filtered_count} molecules")
    logging.info(f"After filtering by fps_pred > {fp_threshold}K: {final_filtered_count} molecules")

    return final_filtered_df

def generate_molecule_image(smiles, output_path):
    """
    Generates and saves a molecule image from a SMILES string.

    Parameters:
    - smiles: str, SMILES string of the molecule.
    - output_path: str, path to save the generated image.
    """
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    opts.bondLineWidth = 2.8

    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = Draw.MolToImage(mol, options=opts)
        img.save(output_path)
    else:
        logging.info(f"Invalid SMILES string: {smiles}")

def ensure_directory_exists(directory_path):
    """
    Ensures the directory exists, if not, creates it.

    Parameters:
    - directory_path: str, path to the directory.
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

# 示例使用
if __name__ == "__main__":
    input_csv_path = 'your_input.csv'
    output_csv_path = 'filtered_molecules.csv'
    mp_threshold = 230  # Replace with your melting point threshold
    bp_threshold = 430  # Replace with your boiling point threshold
    fp_threshold = 360  # Replace with your flash point threshold

    df = read_smiles_from_csv(input_csv_path)
    filtered_df = filter_initial_molecules(df)
    final_filtered_df = filter_molecules(filtered_df, mp_threshold, bp_threshold, fp_threshold)
    final_filtered_df.to_csv(output_csv_path, index=False)
    