import os
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from kpi_features import extract_features

# Disable RDKit warning messages
RDLogger.DisableLog("rdApp.*")

def extract_columns(input_file, output_file):
    """
    Extract specified columns from a CSV file and save to a new CSV file.
    
    Parameters:
    input_file (str): Path to the input CSV file.
    output_file (str): Path to the output CSV file.
    """
    # Read the combined CSV file
    combined_data = pd.read_csv(input_file)
    
    # Extract the specified columns and rename them
    selected_columns = combined_data[['SMILES', 'mp']]
    selected_columns.columns = ['SMILES', 'mp']
    
    # Save the extracted data to a new CSV file
    selected_columns.to_csv(output_file, index=False)
    print(f"Extraction complete and saved to {output_file}")

def standardize_smiles(input_file, output_file):
    """
    Standardize SMILES strings in a CSV file and save to a new CSV file.
    
    Parameters:
    input_file (str): Path to the input CSV file containing SMILES strings.
    output_file (str): Path to the output CSV file for standardized SMILES.
    """
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Create a copy of the DataFrame to save the standardized data
    standardized_df = df.copy()
    
    # Iterate over the SMILES column and standardize each SMILES string
    for index, row in df.iterrows():
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            try:
                standardized_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                standardized_df.at[index, "SMILES"] = standardized_smiles
            except Exception as e:
                print(f"Error processing row {index+1}: {e}")
    
    # Save the standardized data to a new CSV file
    standardized_df.to_csv(output_file, index=False)
    
    # Get the number of rows in the original and new DataFrames
    original_rows = len(df)
    new_rows = len(standardized_df)
    
    print(f"The original table has {original_rows} rows")
    print(f"The new table has {new_rows} rows")

def check_for_duplicates(df, column_name):
    """
    Check for duplicate entries in a specified column of a DataFrame.

    Parameters:
    df (pd.DataFrame): DataFrame to check for duplicates.
    column_name (str): The name of the column to check for duplicates.

    Returns:
    None
    """
    smiles_set = set()
    duplicate_found = False

    # Iterate over the column and check for duplicates
    for index, row in df.iterrows():
        smiles = row[column_name]

        if smiles in smiles_set:
            print(f"Duplicate data found (row {index}): {row}")
            duplicate_found = True
        else:
            smiles_set.add(smiles)

    if not duplicate_found:
        print("Success, no duplicate molecules found")

def filter_and_save_smiles(input_file, output_valid_file, tmp_dir):
    """
    Filters the SMILES in the input CSV file to separate complete and incomplete molecules,
    then saves them to specified output files.

    Parameters:
    - input_file: str, path to the input CSV file.
    - output_valid_file: str, path to the output CSV file for complete molecules.
    - tmp_dir: str, path to the temporary directory for incomplete molecules.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Record the number of rows in the original data
    original_rows = len(df)

    # Filter out SMILES without '.' indicating complete molecules
    valid_smiles_df = df[~df["SMILES"].str.contains('\.')]

    # Filter out SMILES with '.' indicating incomplete molecules
    invalid_smiles_df = df[df["SMILES"].str.contains('\.')]

    # Check and create the temporary directory if it does not exist
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Save complete molecules to a new CSV file
    valid_smiles_df.to_csv(output_valid_file, index=False)

    # Save incomplete molecules to a file in the temporary directory
    output_invalid_file_name = os.path.join(tmp_dir, "1_invalid_molecules.csv")
    invalid_smiles_df.to_csv(output_invalid_file_name, index=False)

    # Get the number of rows in the new dataframe
    new_rows = len(valid_smiles_df)

    # Calculate the number of deleted rows
    deleted_rows = original_rows - new_rows

    # Print the results
    print(f"Original table has {original_rows} rows of data")
    print(f"New table has {new_rows} rows of data")
    print(f"Deleted {deleted_rows} rows containing incomplete molecules with '.'")

def count_elements_in_smiles(input_file, output_file):
    """
    Reads a CSV file containing SMILES strings, counts the occurrence of each element,
    and saves the counts to a new CSV file.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - output_file: str, path to the output CSV file to save the element counts.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Create a dictionary to store the count of each element
    element_counts = {}

    # Iterate through the SMILES column and count elements
    for smiles in df["SMILES"]:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            for atom in mol.GetAtoms():
                element = atom.GetSymbol()
                if element not in element_counts:
                    element_counts[element] = 0
                element_counts[element] += 1

    # Convert the dictionary to a DataFrame for better presentation
    elements_df = pd.DataFrame(list(element_counts.items()), columns=['Element', 'Count'])
    
    # Save the DataFrame to a CSV file
    elements_df.to_csv(output_file, index=False)

    # Print all different elements and their counts
    print("All different elements and their counts:")
    print(elements_df)

def filter_smiles_by_elements(input_file, output_file, allowed_elements):
    """
    Reads a CSV file containing SMILES strings, filters the rows based on allowed elements,
    and saves the filtered data to a new CSV file.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - output_file: str, path to the output CSV file to save the filtered data.
    - allowed_elements: set, set of allowed elements for filtering.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Create an empty DataFrame to store rows that meet the criteria
    selected_rows = pd.DataFrame(columns=df.columns)

    # Iterate through the SMILES column and filter rows based on allowed elements
    for index, row in df.iterrows():
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            if elements.issubset(allowed_elements):
                selected_rows = selected_rows.append(row, ignore_index=True)

    # Save the filtered data to a new CSV file
    selected_rows.to_csv(output_file, index=False)

    # Get the number of rows in the original and new DataFrames
    original_rows = len(df)
    new_rows = len(selected_rows)

    # Print the results
    print(f"Original table has {original_rows} rows of data")
    print(f"New table has {new_rows} rows of data")

def filter_molecules_by_heavy_atoms(input_file, small_output_file, large_output_file):
    """
    Reads a CSV file containing SMILES strings, filters the molecules based on the number of heavy atoms,
    and saves the filtered data to separate CSV files for molecules with <=30 and >30 heavy atoms.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - small_output_file: str, path to the output CSV file for molecules with <=30 heavy atoms.
    - large_output_file: str, path to the output CSV file for molecules with >30 heavy atoms.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Create empty DataFrames to store molecules based on heavy atom count
    selected_rows_small = pd.DataFrame(columns=df.columns)  # Molecules with <=30 heavy atoms
    selected_rows_large = pd.DataFrame(columns=df.columns)  # Molecules with >30 heavy atoms

    # Iterate through the SMILES column and classify molecules based on heavy atom count
    for index, row in df.iterrows():
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            if mol.GetNumHeavyAtoms() <= 30:
                selected_rows_small = selected_rows_small.append(row, ignore_index=True)
            else:
                selected_rows_large = selected_rows_large.append(row, ignore_index=True)

    # Save the filtered data to new CSV files
    selected_rows_small.to_csv(small_output_file, index=False)
    selected_rows_large.to_csv(large_output_file, index=False)

    # Get the number of rows in the original and new DataFrames
    original_rows = len(df)
    small_rows = len(selected_rows_small)
    large_rows = len(selected_rows_large)

    # Print the results
    print(f"Original table has {original_rows} rows of data")
    print(f"Table of molecules with <=30 heavy atoms has {small_rows} rows of data")
    print(f"Table of molecules with >30 heavy atoms has {large_rows} rows of data")

def filter_molecules_by_molecular_weight(input_file, small_mw_output_file, large_mw_output_file, mw_threshold=600):
    """
    Reads a CSV file containing SMILES strings, filters the molecules based on molecular weight,
    and saves the filtered data to separate CSV files for molecules with <=mw_threshold and >mw_threshold molecular weights.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - small_mw_output_file: str, path to the output CSV file for molecules with <=mw_threshold molecular weight.
    - large_mw_output_file: str, path to the output CSV file for molecules with >mw_threshold molecular weight.
    - mw_threshold: float, threshold for molecular weight classification.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Create empty DataFrames to store molecules based on molecular weight
    selected_rows_small_mw = pd.DataFrame(columns=df.columns)  # Molecules with <= mw_threshold molecular weight
    selected_rows_large_mw = pd.DataFrame(columns=df.columns)  # Molecules with > mw_threshold molecular weight

    # Iterate through the SMILES column and classify molecules based on molecular weight
    for index, row in df.iterrows():
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            if mw <= mw_threshold:
                selected_rows_small_mw = selected_rows_small_mw.append(row, ignore_index=True)
            else:
                selected_rows_large_mw = selected_rows_large_mw.append(row, ignore_index=True)

    # Save the filtered data to new CSV files
    selected_rows_small_mw.to_csv(small_mw_output_file, index=False)
    selected_rows_large_mw.to_csv(large_mw_output_file, index=False)

    # Get the number of rows in the original and new DataFrames
    original_rows = len(df)
    small_mw_rows = len(selected_rows_small_mw)
    large_mw_rows = len(selected_rows_large_mw)

    # Print the results
    print(f"Original table has {original_rows} rows of data")
    print(f"Table of molecules with <= {mw_threshold} molecular weight has {small_mw_rows} rows of data")
    print(f"Table of molecules with > {mw_threshold} molecular weight has {large_mw_rows} rows of data")

def filter_molecules_by_property(input_file, property_name, min_value, max_value, satisfied_output_file, not_satisfied_output_file):
    """
    Reads a CSV file containing molecules' data, filters them based on a specified property range,
    and saves the filtered data to separate CSV files for satisfied and not satisfied molecules.

    Parameters:
    - input_file: str, path to the input CSV file containing molecules' data.
    - property_name: str, name of the property to filter by.
    - min_value: float, minimum value for filtering.
    - max_value: float, maximum value for filtering.
    - satisfied_output_file: str, path to the output CSV file for molecules within the specified property range.
    - not_satisfied_output_file: str, path to the output CSV file for molecules outside the specified property range.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Filter molecules that satisfy the property conditions
    df_satisfied = df[(df[property_name] >= min_value) & (df[property_name] <= max_value)]

    # Filter molecules that do not satisfy the property conditions
    df_not_satisfied = df[~((df[property_name] >= min_value) & (df[property_name] <= max_value))]

    # Save the filtered data to new CSV files
    df_satisfied.to_csv(satisfied_output_file, index=False)
    df_not_satisfied.to_csv(not_satisfied_output_file, index=False)

    # Output the number of molecules in each file
    num_satisfied = df_satisfied.shape[0]
    num_not_satisfied = df_not_satisfied.shape[0]

    print(f"Number of molecules satisfying the condition: {num_satisfied}")
    print(f"Number of molecules not satisfying the condition: {num_not_satisfied}")

def process_smiles(input_file, success_output_file, error_output_file):
    """
    Reads a CSV file containing SMILES strings, processes each SMILES using the extract_features function,
    and saves the successful and failed results to separate CSV files.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - success_output_file: str, path to the output CSV file for successfully processed SMILES.
    - error_output_file: str, path to the output CSV file for SMILES that caused errors during processing.
    """
    # Read the original CSV file
    df = pd.read_csv(input_file)

    # Create DataFrames to store successful and failed data
    success_df = pd.DataFrame(columns=df.columns)
    error_df = pd.DataFrame(columns=df.columns)

    # Iterate through each row, extract and process the SMILES column
    for index, row in df.iterrows():
        try:
            # Extract the SMILES column information
            smi = row['SMILES']

            # Call the extract_features function to process the SMILES
            extract_features(smi)

            # If no error occurs, save the row to the success DataFrame
            success_df = success_df.append(row, ignore_index=True)
        except Exception as e:
            # If an error occurs, save the row to the error DataFrame
            error_df = error_df.append(row, ignore_index=True)

    # Save the successful and failed data to new CSV files
    success_df.to_csv(success_output_file, index=False)
    error_df.to_csv(error_output_file, index=False)

    # Print the number of rows in the original data and the two CSV files
    print("Number of rows in original data:", len(df))
    print("Number of rows in successful data:", len(success_df))
    print("Number of rows in failed data:", len(error_df))
    